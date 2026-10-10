/* Crackle test harness: evolve a single gas parcel for nsteps system timesteps.
 *
 * usage: crackle_harness use_grackle nH T0 Z uvb isrf nsteps dt_myr nsub ssm h2ss dustevol
 *   use_grackle  2 = crackle, 1 = original Fortran grackle
 *   nH           hydrogen number density [cm^-3]
 *   T0           initial temperature [K]
 *   Z            metallicity [solar]
 *   uvb          UVbackground flag
 *   isrf         ISRF in Habing units; < 0 marks a non-ISM particle (crackle ism_flag=0)
 *   nsteps       number of system timesteps
 *   dt_myr       system timestep [Myr]
 *   nsub         split each system timestep into nsub calls (nsub>>1 gives a convergence reference)
 *   ssm          self_shielding_method
 *   h2ss         H2_self_shielding
 *   dustevol     use_dust_evol
 *
 * environment:
 *   CRACKLE_DATA_DIR  directory containing the CloudyData_*.h5 files (default: ".")
 *   ACC               chemistry->accuracy (default 0.1)
 *   RTH2              if set, turns on radiative transfer with this H2 dissociation rate [1/s]
 *
 * output, one line per system timestep: t[Myr] T[K] xHI xe xH2 dust2gas
 * (xHI, xe, xH2 are mass fractions relative to total hydrogen)
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <grackle.h>

#define ALLOC(f) f = calloc(1, sizeof(gr_float))

int main(int argc, char **argv)
{
  if (argc < 13) {
    fprintf(stderr, "usage: %s use_grackle nH T0 Z uvb isrf nsteps dt_myr nsub ssm h2ss dustevol\n", argv[0]);
    return 1;
  }
  const int use_grackle = atoi(argv[1]);
  const double nH = atof(argv[2]), T0 = atof(argv[3]), Z = atof(argv[4]);
  const int uvb = atoi(argv[5]);
  const double isrf = atof(argv[6]);
  const int nsteps = atoi(argv[7]);
  const double dt_myr = atof(argv[8]);
  const int nsub = atoi(argv[9]);
  const int ssm = atoi(argv[10]), h2ss = atoi(argv[11]), dustevol = atoi(argv[12]);
  const double accuracy = getenv("ACC") ? atof(getenv("ACC")) : 0.1;
  const double rth2 = getenv("RTH2") ? atof(getenv("RTH2")) : -1.;  /* < 0 => RT off */
  const char *datadir = getenv("CRACKLE_DATA_DIR") ? getenv("CRACKLE_DATA_DIR") : ".";

  /* Physical units at z=0; density unit ~ m_H so code density ~ number density */
  grackle_verbose = 0;
  code_units u;
  u.comoving_coordinates = 0;
  u.density_units = 1.67e-24;
  u.length_units = 3.086e21;
  u.time_units = 3.15e13;
  u.a_units = 1.0;
  u.a_value = 1.;
  set_velocity_units(&u);
  const double tunits = get_temperature_units(&u);
  /* note: u.temperature_units deliberately not set, to test that calculate_temperature works before solve_chemistry */

  chemistry_data *cd = malloc(sizeof(chemistry_data));
  set_default_chemistry_parameters(cd);
  grackle_data->use_grackle = use_grackle;
  grackle_data->with_radiative_cooling = 1;
  grackle_data->primordial_chemistry = 2;
  grackle_data->dust_chemistry = 1;
  grackle_data->use_dust_evol = dustevol;
  grackle_data->use_dust_density_field = dustevol;
  grackle_data->metal_cooling = 1;
  grackle_data->UVbackground = uvb;
  grackle_data->self_shielding_method = ssm;
  grackle_data->H2_self_shielding = h2ss;
  grackle_data->accuracy = accuracy;
  grackle_data->max_iterations = 10000;
  grackle_data->use_isrf_field = 1;
  grackle_data->use_radiative_transfer = rth2 >= 0.;
  grackle_data->sne_coeff = 100;
  grackle_data->dust_growth_densref = 2.3e-27;
  static char datafile[4096];
  snprintf(datafile, sizeof(datafile), "%s/%s", datadir, ssm > 0 ? "CloudyData_UVB=HM2012_shielded.h5" : "CloudyData_UVB=HM2012.h5");
  grackle_data->grackle_data_file = datafile;
  const double sa[10] = {0.2485, 2.38e-3, 0.70e-3, 5.79e-3, 1.26e-3, 7.14e-4, 6.71e-4, 3.12e-4, 0.65e-4, 1.31e-3};
  for (int i = 0; i < 10; i++) grackle_data->SolarAbundances[i] = sa[i];
  if (initialize_chemistry_data(&u) == 0) {
    fprintf(stderr, "initialize_chemistry_data failed (data file %s)\n", datafile);
    return 1;
  }

  /* Single-element fields (crackle works on one particle per call) */
  grackle_field_data f;
  f.grid_rank = 3;
  f.grid_dimension = malloc(3*sizeof(int));
  f.grid_start = malloc(3*sizeof(int));
  f.grid_end = malloc(3*sizeof(int));
  for (int i = 0; i < 3; i++) { f.grid_dimension[i] = 1; f.grid_start[i] = 0; f.grid_end[i] = 0; }
  f.grid_dx = 0.;
  ALLOC(f.density); ALLOC(f.internal_energy); ALLOC(f.x_velocity); ALLOC(f.y_velocity); ALLOC(f.z_velocity);
  ALLOC(f.HI_density); ALLOC(f.HII_density); ALLOC(f.HM_density);
  ALLOC(f.HeI_density); ALLOC(f.HeII_density); ALLOC(f.HeIII_density);
  ALLOC(f.H2I_density); ALLOC(f.H2II_density); ALLOC(f.DI_density); ALLOC(f.DII_density); ALLOC(f.HDI_density);
  ALLOC(f.e_density); ALLOC(f.metal_density); ALLOC(f.volumetric_heating_rate); ALLOC(f.specific_heating_rate);
  f.temperature_floor = NULL;
  ALLOC(f.RT_HI_ionization_rate); ALLOC(f.RT_HeI_ionization_rate); ALLOC(f.RT_HeII_ionization_rate);
  ALLOC(f.RT_H2_dissociation_rate); ALLOC(f.RT_heating_rate);
  ALLOC(f.H2_self_shielding_length); ALLOC(f.H2_custom_shielding_factor); ALLOC(f.isrf_habing);
  ALLOC(f.dust_density); ALLOC(f.SNe_ThisTimeStep);
  gr_float **gas_met[10] = {&f.He_gas_metalDensity, &f.C_gas_metalDensity, &f.N_gas_metalDensity, &f.O_gas_metalDensity, &f.Ne_gas_metalDensity,
                            &f.Mg_gas_metalDensity, &f.Si_gas_metalDensity, &f.S_gas_metalDensity, &f.Ca_gas_metalDensity, &f.Fe_gas_metalDensity};
  gr_float **dust_met[10] = {&f.He_dust_metalDensity, &f.C_dust_metalDensity, &f.N_dust_metalDensity, &f.O_dust_metalDensity, &f.Ne_dust_metalDensity,
                             &f.Mg_dust_metalDensity, &f.Si_dust_metalDensity, &f.S_dust_metalDensity, &f.Ca_dust_metalDensity, &f.Fe_dust_metalDensity};
  for (int i = 0; i < 10; i++) { ALLOC(*gas_met[i]); ALLOC(*dust_met[i]); }

  /* Initial state: ionization set roughly by T0 so we don't start wildly out of equilibrium */
  const double X = grackle_data->HydrogenFractionByMass;
  const double rho = nH / X;
  const double tiny = 1.e-20;
  const double xion = T0 > 2.e4 ? 0.999 : (T0 > 8.e3 ? 0.5 : 1.e-4);
  f.density[0] = rho;
  f.HI_density[0] = X*rho*(1.-xion);
  f.HII_density[0] = X*rho*xion;
  f.HM_density[0] = tiny*rho;
  f.H2I_density[0] = (T0 < 1.e3 ? 0.1 : 1.e-6)*X*rho;
  f.HI_density[0] -= f.H2I_density[0];
  f.H2II_density[0] = tiny*rho;
  f.HeI_density[0] = (1.-X)*rho*(1.-xion);
  f.HeII_density[0] = (1.-X)*rho*xion*0.5;
  f.HeIII_density[0] = (1.-X)*rho*xion*0.5;
  f.DI_density[0] = 2.*3.4e-5*rho; f.DII_density[0] = tiny*rho; f.HDI_density[0] = tiny*rho;
  f.e_density[0] = f.HII_density[0] + 0.25*f.HeII_density[0] + 0.5*f.HeIII_density[0];

  /* Metals; with dust evolution, 30% of non-He metals start in dust */
  const double d2m = dustevol ? 0.3 : 0.;
  f.metal_density[0] = rho * grackle_data->SolarMetalFractionByMass * Z * (1.-d2m);
  f.dust_density[0] = 0.;
  for (int i = 0; i < 10; i++) {
    const double tot = rho*sa[i]*Z;
    const double dfrac = (i == 0) ? 0. : d2m;
    (*gas_met[i])[0] = tot*(1.-dfrac);
    (*dust_met[i])[0] = tot*dfrac;
    f.dust_density[0] += tot*dfrac;
  }
  f.SNe_ThisTimeStep[0] = 0.;
  f.isrf_habing[0] = isrf;
  if (rth2 >= 0.) f.RT_H2_dissociation_rate[0] = rth2 * u.time_units;

  /* Internal energy for T0 */
  const double ntot = f.HI_density[0] + f.HII_density[0] + 0.25*(f.HeI_density[0] + f.HeII_density[0] + f.HeIII_density[0])
                      + 0.5*f.H2I_density[0] + f.e_density[0];
  f.internal_energy[0] = T0 / (rho / ntot * (grackle_data->Gamma - 1.) * tunits);

  const double dt = dt_myr * 3.15e13 / u.time_units;
  gr_float T[1];
  for (int s = 0; s <= nsteps; s++) {
    if (s > 0) {
      for (int k = 0; k < nsub; k++) {
        if (solve_chemistry(&u, &f, dt/nsub) == 0) { fprintf(stderr, "solve_chemistry failed\n"); return 1; }
      }
    }
    calculate_temperature(&u, &f, T);
    printf("%10.4f %12.5e %10.4e %10.4e %10.4e %10.4e\n", s*dt_myr, T[0], f.HI_density[0]/(X*rho),
           f.e_density[0]/(X*rho), f.H2I_density[0]/(X*rho), f.dust_density[0]/rho);
  }
  return 0;
}
