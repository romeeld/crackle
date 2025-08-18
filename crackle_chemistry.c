
/* Main driver for Crackle chemistry solver, which is a fully C version of Grackle .
 * The code explicitly evolves the internal energy, electron density, and HI density.
 * Other quantities are evolved via a predictor-corrector scheme, if explicit evolution is not accurate enough.
 *
 %g %g * Romeel Dave, Feb 2024
 * */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "phys_constants.h"
/* Grackle includes */
//#include "grackle.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
/* Crackle-specific includes */
#include "crackle.h"
#include "crackle_aux.h"
#include "crackle_interpolation.h"
#include "crackle_dust.h"
//#include <omp.h>

#define MINFRAC 1.e-3
#define CONVERGENCE 1.e-2
#define HEATLIM 10.  // max heating factor allowed in a single system timestep


int crackle_solve_chemistry(grackle_field_data *p, chemistry_data *chemistry, chemistry_data_storage grackle_rates, photo_rate_storage my_uvb_rates, code_units *units, double dt) 
{
	int iter=0;
	double dtit=0., dtcool=0., dtsuppress=1.; 
	grackle_part_data gp, gp_old;
	chemistry_rate_storage my_rates;  // interpolated rates for this field
	interp_struct interpolation;
	crackle_units cunits;

	//if (chemistry->UVbackground == 0) {
	/* Set up uv rates to be 0. */
	//	my_uvb_rates.k24 = my_uvb_rates.k25 = my_uvb_rates.k26 =
        //	my_uvb_rates.k27 = my_uvb_rates.k28 = my_uvb_rates.k29 =
        //	my_uvb_rates.k30 = my_uvb_rates.k31 = my_uvb_rates.piHI =
        //	my_uvb_rates.piHeI = my_uvb_rates.piHeII = my_uvb_rates.crsHI =
        //	my_uvb_rates.crsHeI = my_uvb_rates.crsHeII =
        //	my_uvb_rates.comp_xray = my_uvb_rates.temp_xray = 0.;

		//my_uvb_rates.k26 = 0.0;
	//	my_rates.k24 = 0.0;
	//	my_rates.k25 = 0.0;
	//	my_rates.k26 = 0.0;
	//	my_rates.k27 = 0.0;
	//	my_rates.k28 = 0.0;
	//	my_rates.k29 = 0.0;
	//	my_rates.k30 = 0.0;
	//	my_rates.k31 = 0.0;

	//	my_rates.piHI      = 0.0;
	//	my_rates.piHeI     = 0.0;
	//	my_rates.piHeII    = 0.0;
	//	my_rates.crsHI     = 0.0;
	//	my_rates.crsHeI    = 0.0;
	//	my_rates.crsHeII   = 0.0;
	//	my_rates.comp_xray = 0.0;
	//}

	/* Copy to grackle_part_data */
	copy_grackle_fields_to_part(p, &gp, chemistry);
	
	/* Checking Nan for gp field data*/
	//if (chemistry->use_radiative_transfer){
	//	if (has_nan_in_gp(&gp)) {
    	//		printf("One or more NaNs found in gp structure before cooling calculation.\n");
	//	}
	//}
	
	if (chemistry->use_radiative_transfer) {
    		FILE *logfile = fopen("gp_nan_log.txt", "a");
    		if (!logfile) perror("Failed to open log file");

    		int fixed = fix_nan_in_gp(&gp, logfile);
    		if (fixed > 0) {
        		printf("%d NaN values replaced with 0 in gp structure.\n", fixed);
    		}

    		if (logfile) fclose(logfile);
	}

	/* Set up various unit conversions etc */
	set_crackle_units(units, grackle_rates, chemistry->Gamma, &cunits);

	/* Compute some basic properties */
	set_rhot(&gp, units, chemistry);

	/* initialize inteprolation of chemistry rate tables for this particle */
	init_temperature_interpolation(&gp, chemistry, &interpolation, cunits, grackle_rates);

	/* ISM flag: see whether we should be evolving dust and H2 */
	int ism_flag = (gp.isrf_habing >= 0.);
	if (gp.isrf_habing < 0.) gp.isrf_habing = 0.;

	//printf("=== Crackle RT ionization & heating rates ===\n");
        //printf("Crackle RT_heating_rate         = %e\n", p->RT_heating_rate[0]);
        //printf("Crackle RT_HI_ionization_rate   = %e\n", p->RT_HI_ionization_rate[0]);
        //printf("Crackle RT_HeI_ionization_rate  = %e\n", p->RT_HeI_ionization_rate[0]);
        //printf("Crackle RT_HeII_ionization_rate = %e\n", p->RT_HeII_ionization_rate[0]);
        //printf("Crackle RT_H2_dissociation_rate = %e\n", p->RT_H2_dissociation_rate[0]);

	while (dtcool < dt) {
	    //if (gp.density > 2.e9) gp.verbose=1;
	    /* Retain previous iteration particle info */
	    memmove(&gp_old, &gp, sizeof(gp)); 
	    /* Set up cooling/heating rates interpolation */
	    setup_temperature_interpolation(0.5*(gp.tgas+gp_old.tgas), chemistry, &interpolation);
	    /* Get interpolated chemistry rates for this particle */
	    lookup_chemistry_coeffs(chemistry->primordial_chemistry, grackle_rates, &my_rates, &interpolation);  
	    /* Compute rate of change of thermal energy */
	    compute_edot(&gp, chemistry, grackle_rates, &my_rates, my_uvb_rates, &interpolation, units, cunits);  
	    /* If we've reached temp floor and are still cooling, then we're done */
	    if (apply_temperature_bounds(&gp, chemistry, gp.temperature_floor, HEATLIM * gp.tgas)) break;

	    /* Check if we can (and if it's worthwhile to) use a predictor-corrector for the rest of the timestep */
	    if (!ism_flag && fabs(gp.edot_ext) > EDOT_EXT_FACTOR * fabs(gp.edot-gp.edot_ext) && dtit < 0.25 * (dt-dtcool)) {
		/* Predictor step */
		dtit = dt - dtcool;  // do entire remainder of the timestep at once
	        evolve_hydrogen(&gp, &gp, chemistry, my_rates, dtit);  
	        evolve_helium(&gp, &gp, chemistry, my_rates, dtit); 
	        if (chemistry->primordial_chemistry >= 2) evolve_H2(&gp, ism_flag, chemistry, my_rates, cunits, dtit);  // sets H2=0
	        evolve_elements(&gp, &gp_old, chemistry);
	        evolve_internal_energy(&gp, chemistry, dtit);
	        set_rhot(&gp, units, chemistry);
	        memmove(&gp_old, &gp, sizeof(gp)); 
	        compute_edot(&gp, chemistry, grackle_rates, &my_rates, my_uvb_rates, &interpolation, units, cunits);  
	        evolve_hydrogen(&gp, &gp, chemistry, my_rates, dtit);  
	        evolve_helium(&gp, &gp, chemistry, my_rates, dtit); 
	        evolve_elements(&gp, &gp_old, chemistry);
	        evolve_internal_energy(&gp, chemistry, dtit);
	        evolve_pred_corr(&gp, &gp_old, chemistry);
	        if (chemistry->use_dust_evol) evolve_dust(&gp, chemistry, units, ism_flag, dtit); 
		dtcool = dt + tiny;
	        apply_temperature_bounds(&gp, chemistry, gp.temperature_floor, HEATLIM * gp.tgas);
		break;
	    }

	    /* Compute rate of change of electron density, used to set the timestep */
	    gp.HIdot = compute_HIdot(chemistry->primordial_chemistry, gp, chemistry, my_rates, units); 
	    gp.dedot = compute_dedot(chemistry->primordial_chemistry, gp, chemistry, my_rates, units); 
	    /* Set the timestep for this iteration (lowering rates if they are too large) */
	    dtit = compute_iteration_dt(&gp, &gp_old, chemistry, dt, dtcool, &dtsuppress);

	    /* Evolve all quantities */
	    evolve_hydrogen(&gp, &gp, chemistry, my_rates, dtit);  
	    evolve_helium(&gp, &gp, chemistry, my_rates, dtit); 
	    if (chemistry->primordial_chemistry >= 2) evolve_H2(&gp, ism_flag, chemistry, my_rates, cunits, dtit);  
	    if (chemistry->use_dust_evol) evolve_dust(&gp, chemistry, units, ism_flag, dtit); 
	    evolve_elements(&gp, &gp_old, chemistry);
	    /* Advance thermal energy over full step using midpoint edot */
	    evolve_internal_energy(&gp, chemistry, dtit);
	    set_rhot(&gp, units, chemistry);
	    dtcool += dtit;
	    iter ++;

	    if (gp.verbose) printf("iter: i=%d dt=%g nh=%g e=%g de=%g HI=%g HII=%g H2I=%g HeI=%g HeII=%g dust=%g T=%g Td=%g\n", iter, dtit, gp.rhoH * units->density_units / mh, gp.internal_energy, gp.e_density, gp.HI_density/gp.density, gp.HII_density/gp.density, gp.H2I_density/gp.density, gp.HeI_density/gp.density, gp.HeII_density/gp.density, gp.dust_density/gp.density, gp.tgas, gp.tdust);
	    if (gp.verbose) printf("rates: i=%d fdt=%g edot=%g dedot=%g HIdot=%g T=%g n=%g\n",iter, dtit/dt, gp.edot, gp.dedot, gp.HIdot, gp.tgas, gp.density*units->density_units/mh );
	    
	    /* Check for convergence or too many iterations */
	    //if (fabs(gp.edot * dtit) < CONVERGENCE * gp.internal_energy * gp.density && fabs(gp.HIdot * dtit) < CONVERGENCE * fmax(gp.HI_density, MINFRAC) && fabs(gp.dedot * dtit) < CONVERGENCE * fmax(gp.e_density, MINFRAC)) break;
	    if (fabs(gp_old.internal_energy - gp.internal_energy) < CONVERGENCE * gp_old.internal_energy 
		&& fabs(gp_old.HI_density - gp.HI_density) < CONVERGENCE * fmax(gp_old.HI_density, MINFRAC*gp_old.density)
		&& fabs(gp_old.e_density - gp.e_density) < CONVERGENCE * fmax(gp_old.e_density, MINFRAC*gp_old.density) ) break;
	    if (iter > chemistry->max_iterations) break;
	    gp.verbose = 0;
	}

	/* If there is time left over, need to form/destroy dust over the full timestep */
	if (dtcool < dt && chemistry->primordial_chemistry >= 2) {
	    evolve_H2(&gp, ism_flag, chemistry, my_rates, cunits, dt - dtcool); 
	    evolve_elements(&gp, &gp_old, chemistry);
	}
	if (dtcool < dt && chemistry->use_dust_evol) evolve_dust(&gp, chemistry, units, ism_flag, dt - dtcool); 

	/* Copy from grackle_part_data */
	copy_grackle_fields_from_part(p, &gp, chemistry);

	return 1;
}

double compute_iteration_dt(grackle_part_data *gp, grackle_part_data *gp_old, chemistry_data *chemistry, double dt, double dtcool, double *dtsuppress)
{
	/* Calculate timesteps for each tracked quantity */
	const double dt_e = fabs( gp->internal_energy * gp->density / (gp->edot + tiny) );
	const double dt_de = fabs( fmax(gp->e_density, MINFRAC*gp->density) / (gp->dedot + tiny) );
	const double dt_HI = fabs( fmax(gp->HI_density, MINFRAC*gp->density) / (gp->HIdot + tiny) );

	/* If thermal equilbirium is passed, attenuate the timestep to help it converge */
	if (gp->edot * gp_old->edot < 0.f) *dtsuppress *= 0.5;
	else *dtsuppress *= (1. + 0.5*chemistry->accuracy);  // asymmetric wrt the reduction to avoid oscillations

	/* Compute the timestep for this iteration */
	double dtit = chemistry->accuracy * *dtsuppress * fmin(fmin(dt_de, dt_HI), dt_e);

	/* Limit dtit if we will reach end of timestep this iteration */
	if (dtit > dt - dtcool) {
	    dtit = dt - dtcool + tiny;
	}

	return dtit;
}

void evolve_internal_energy(grackle_part_data *gp, chemistry_data *chemistry, double dtit) 
{
	const double u_prev = gp->internal_energy;
	gp->internal_energy += gp->edot / gp->density * dtit;
	/* Limits -- don't let u change by more than accuracy level in a single iteration */
	if (gp->internal_energy < (1.-chemistry->accuracy) * u_prev) gp->internal_energy = (1.-chemistry->accuracy) * u_prev;
	if (gp->internal_energy > (1.+chemistry->accuracy) * u_prev) gp->internal_energy = (1.+chemistry->accuracy) * u_prev;
	if (gp->internal_energy < gp->u_cmb) gp->internal_energy = gp->u_cmb;

	return;
}

void compute_edot(grackle_part_data *gp, chemistry_data *chemistry, chemistry_data_storage grackle_rates, chemistry_rate_storage *my_rates, photo_rate_storage my_uvb_rates, interp_struct *interpolation, code_units *units, crackle_units cunits)
{
	double edot_prim = 0., edot_h2 = 0., edot_gasgr = 0., edot_uvb = 0., edot_pe = 0., edot_edust = 0., edot_comp = 0., edot_rt = 0., edot_h2heat = 0., edot_ext = 0., edot_metal = 0.;

	/* Interpolate primordial excitation, ionization, and recombination rates for given T */
	primordial_cooling_rates(interpolation, grackle_rates, my_rates);

	gp->edot = 0.;  // initialize for this step

	/* Primordial cooling */
	edot_prim += -my_rates->ceHI*gp->HI_density*gp->e_density 
		- 0.25*my_rates->ceHeI*gp->HeII_density*gp->e_density*gp->e_density*cunits.dom 
		- 0.25*my_rates->ceHeII*gp->HeII_density*gp->e_density;
	edot_prim += -my_rates->ciHI*gp->HI_density*gp->e_density
		- 0.25*my_rates->ciHeI*gp->HeI_density*gp->e_density
		- 0.25*my_rates->ciHeII*gp->HeII_density*gp->e_density
		- 0.25*my_rates->ciHeIS*gp->HeII_density*cunits.dom*gp->e_density*gp->e_density;
	edot_prim += -my_rates->reHII* gp->HII_density*gp->e_density
		- 0.25*my_rates->reHeII1*gp->HeII_density*gp->e_density
		- 0.25*my_rates->reHeII2*gp->HeII_density*gp->e_density
		- 0.25*my_rates->reHeIII*gp->HeIII_density*gp->e_density;
	edot_prim -= my_rates->brem*(gp->HII_density+0.25*gp->HeII_density+gp->HeIII_density)*gp->e_density;
	gp->edot += edot_prim;

	/* Add H2 cooling (Glover & Abel 2008); fudge at extreme density NOT included */
	if (chemistry->primordial_chemistry >= 2) {
	    H2_cooling_rates(interpolation, grackle_rates, my_rates);
	    const double galdl = my_rates->GAHI * gp->HI_density + 0.5 * my_rates->GAH2 * gp->H2I_density + 0.25 * my_rates->GAHe * gp->HeI_density + my_rates->GAHp * gp->HII_density + my_rates->GAel * gp->e_density;
	    const double gphdl1 = my_rates->H2LTE * cunits.dom_inv;
	    edot_h2 -= gp->H2I_density * my_rates->H2LTE / ((1.+gphdl1/galdl) * 2.*cunits.dom);
	    gp->edot += edot_h2;

	    /* Dust-related cooling/heating processes */
	    if (chemistry->use_dust_evol || chemistry->dust_chemistry) {
	        /* compute dust temperature */
	        double gasgr = interpolate_rates(grackle_rates.gas_grain, interpolation->binfrac, interpolation->index);
	        double gasgr_tdust = chemistry->local_dust_to_gas_ratio * gasgr * cunits.coolunit / mh;
	        /* Get dust rates */
	        //gp->tdust = calculate_dust_temp(gp->tgas, gp->nH, gasgr_tdust, grackle_rates.gamma_isrf * gp->isrf_habing, cunits.compton2, gp->tdust); // calc_tdust_1d_g()
	        calculate_tdust_bisect(gp, gasgr_tdust, grackle_rates.gamma_isrf * gp->isrf_habing, cunits.compton2);
		//if(gp->tdust < 2*cunits.compton2) printf("DUST td=%g gasgr=%g tcmb=%g isrf=%g\n",gp->tdust, gasgr_tdust, cunits.compton2, gp->isrf_habing);
	        /* Gas-dust grain heat transfer rate */
	        dust_species_rates(gp->tdust, gp->dust2gas, chemistry, grackle_rates, my_rates, interpolation);
	        if (gp->dust2gas > tiny) {
	            edot_gasgr -= gasgr * (gp->tgas - gp->tdust) * gp->dust2gas * gp->rhoH * gp->rhoH;
	            gp->edot += edot_gasgr;
	        }
	    }
	}

	/* UVB heating and shielding */
	if (chemistry->UVbackground && gp->grid_end >= 0) {
	    compute_self_shielded_rates(gp, chemistry, my_rates, my_uvb_rates, cunits);
	    edot_uvb += cunits.dom_inv * (my_uvb_rates.piHI * gp->fSShHI * gp->HI_density + 0.25 * (my_uvb_rates.piHeI * gp->fSShHeI * gp->HeI_density + my_uvb_rates.piHeII * gp->fSShHeII * gp->HeII_density));
	    gp->edot += edot_uvb;
	}

	/* Photo-electric heating by UV-irradiated dust */
	if (chemistry->use_isrf_field && gp->isrf_habing > tiny && gp->dust2gas > tiny) {
	    double gammah_eff = my_rates->gammah * 0.05 * gp->isrf_habing;  // Wolfire+1995 eq 1
	    if (gp->tgas > 2.e4) gammah_eff = 0.;
	    edot_pe += gammah_eff * gp->rhoH * cunits.dom_inv * gp->dust2gas / chemistry->local_dust_to_gas_ratio;
	    gp->edot += edot_pe;
	    //printf("pe %g %g %g\n",gp->edot,gp->dust2gas,gammah_eff);
	}

	/* Electron recombination onto dust grains */
	if ((chemistry->use_dust_evol || chemistry->dust_chemistry) && gp->metallicity > tiny && gp->e_density > tiny && gp->isrf_habing > tiny) {
	    const double regr = interpolate_rates(grackle_rates.regr, interpolation->binfrac, interpolation->index);
	    const double grbeta = 0.74 * pow(gp->tgas, -0.068);
	    edot_edust -= regr * pow(gp->isrf_habing * cunits.dom_inv / gp->e_density, grbeta) * gp->e_density * gp->rhoH * gp->dust2gas / chemistry->local_dust_to_gas_ratio;
	    gp->edot += edot_edust;
	}

	/* CMB Compton cooling and X-ray compton heating */
	edot_comp -= cunits.compton1 * (gp->tgas-cunits.compton2) * gp->e_density * cunits.dom_inv;
	edot_comp -= my_uvb_rates.comp_xray * (gp->tgas-my_uvb_rates.temp_xray) * gp->e_density * cunits.dom_inv;
	gp->edot += edot_comp;

	/* Photoheating from radiative transfer */
	if (chemistry->use_radiative_transfer && gp->RT_heating_rate > tiny) {
	    edot_rt += gp->RT_heating_rate * gp->HI_density * cunits.dom_inv * cunits.coolunit_inv;
	    
	    // === NaN Check ===
    	if (isnan(edot_rt) || isnan(gp->RT_heating_rate) || isnan(gp->HI_density) ||
        	isnan(cunits.dom_inv) || isnan(cunits.coolunit_inv)) {
        	printf("NaN detected!\n");
        	printf("RT_heating_rate: %g\n", gp->RT_heating_rate);
        	printf("HI_density: %g\n", gp->HI_density);
        	printf("dom_inv: %g\n", cunits.dom_inv);
        	printf("coolunit_inv: %g\n", cunits.coolunit_inv);
        	printf("edot_rt: %g\n", edot_rt);
        	exit(1);  // Optional: stop the program immediately
    	}

	    gp->edot += edot_rt;
	    gp->edot_ext += edot_rt;
	}

	/* H2 heating */
	if (chemistry->primordial_chemistry >= 2) {
	    double h2heatfac = 1. + my_rates->n_cr_n / (cunits.dom * my_rates->n_cr_d1 * gp->HI_density + 0.5 * my_rates->n_cr_d2 * gp->H2I_density);
	    h2heatfac = 1./h2heatfac;
	    double H2delta = 4.48 * gp->HI_density * (my_rates->k22 * gp->HI_density * gp->HI_density - 0.5 * my_rates->k13 * gp->H2I_density);
	    if (H2delta > tiny) H2delta *= h2heatfac;
	    if (chemistry->use_dust_evol && gp->dust2gas > tiny) {
		H2delta += my_rates->h2dust * gp->HI_density * gp->rhoH * (0.2 + 4.2 * h2heatfac);
	    }
	    edot_h2heat += cunits.chunit * H2delta;
	    gp->edot += edot_h2heat;
	}

	/* Cloudy-based metal cooling/heating */
	if (chemistry->metal_cooling && gp->metallicity > tiny) {
	    edot_metal += gp->rhoH * gp->rhoH * gp->metallicity * cloudy_metal_cooling(gp->logtem, chemistry, interpolation, grackle_rates, cunits);  // cool1d_cloudy_g()
	    gp->edot += edot_metal;
	}

	/* User-specified external heating terms */
	if (chemistry->use_volumetric_heating_rate && gp->volumetric_heating_rate > tiny) {
	    edot_ext += gp->volumetric_heating_rate * cunits.dom_inv * cunits.dom_inv * cunits.coolunit_inv;
	}
	if (chemistry->use_specific_heating_rate && gp->specific_heating_rate > tiny) {
	    edot_ext += gp->specific_heating_rate * gp->density * mh * cunits.dom_inv * cunits.coolunit_inv;
	}
	gp->edot += edot_ext;
	gp->edot_ext += edot_ext;

	//gp->edot = edot_prim + edot/_h2 + edot_gasgr + edot_uvb + edot_pe + edot_edust + edot_comp + edot_rt + edot_h2heat + edot_ext + edot_metal;
	if (gp->verbose) {
	    printf("edot: %g pr=%g h2=%g gr=%g uvb=%g pe=%g ed=%g co=%g rt=%g h2h=%g ext=%g met=%g\n",gp->edot, edot_prim , edot_h2 , edot_gasgr , edot_uvb, edot_pe , edot_edust , edot_comp , edot_rt , edot_h2heat , edot_ext , edot_metal); 
	    fflush(stdout);
	}

	assert(gp->edot==gp->edot && "CRACKLE ERROR: gp->edot is nan.");  // check for NaN

	return;
}

double compute_dedot(int chemistry_flag, grackle_part_data gp, chemistry_data *chemistry, chemistry_rate_storage my_rates, code_units *units) 
{ 
	double dedot=0.;

	/* Recombination */
	dedot = my_rates.k1 * gp.HI_density * gp.e_density 
		+ 0.25 * my_rates.k3 * gp.HeI_density * gp.e_density 
		+ 0.25 * my_rates.k5 * gp.HeII_density * gp.e_density
		- my_rates.k2 * gp.HII_density * gp.e_density
		- 0.25 * my_rates.k4 * gp.HeII_density * gp.e_density
		- 0.25 * my_rates.k6 * gp.HeIII_density * gp.e_density
		+ my_rates.k57 * gp.HI_density * gp.HI_density 
		+ 0.25 * my_rates.k58 * gp.HI_density * gp.HeI_density;

	if (chemistry->UVbackground > 0) {
		+ my_rates.k24shield * gp.HI_density 
		+ 0.25 * my_rates.k25shield * gp.HeII_density 
		+ 0.25 * my_rates.k26shield * gp.HeI_density;
	}
	/* RT photoionization */
	if (chemistry->use_radiative_transfer && gp.RT_heating_rate > tiny) {
	    dedot += gp.RT_HI_ionization_rate * gp.HI_density + 0.25 * (gp.RT_HeI_ionization_rate * gp.HeI_density + gp.RT_HeII_ionization_rate * gp.HeII_density);
	}

	/* Molecular */
	if (chemistry_flag > 1) {
	    dedot += my_rates.k8 * gp.HI_density * gp.HM_density
		+ my_rates.k15 * gp.HI_density * gp.HM_density
		+ my_rates.k17 * gp.HII_density * gp.HM_density
		+ my_rates.k14 * gp.e_density * gp.HM_density
	    	- my_rates.k7 * gp.HI_density * gp.e_density 
		- 0.5 * my_rates.k18 * gp.H2II_density * gp.e_density;
	}

	return dedot;
}

double compute_HIdot(int chemistry_flag, grackle_part_data gp, chemistry_data *chemistry, chemistry_rate_storage my_rates, code_units *units) 
{ 
	double HIdot=0.;

	/* Recombination */
	HIdot = -my_rates.k1* gp.HI_density * gp.e_density 
		+ my_rates.k2* gp.HII_density * gp.e_density 
	/* Collisional ionization */
		- my_rates.k57 * gp.HI_density * gp.HI_density 
		- 0.25 * my_rates.k58 * gp.HI_density * gp.HeI_density 
	/* UVB photo-ionization (no self-shielding -- assumes self-shielding built into ionization table if desired */
		- my_rates.k24shield * gp.HI_density;
	/* RT photoionization */
	if (chemistry->use_radiative_transfer) {
	    HIdot -= gp.RT_HI_ionization_rate * gp.HI_density;
	}

	/* Molecular */
	if (chemistry_flag >= 2) {
	    HIdot += -my_rates.k7 * gp.e_density * gp.HI_density 
		    - my_rates.k8 * gp.HM_density * gp.HI_density 
		    - my_rates.k9 * gp.HII_density * gp.HI_density 
		    - 0.5 *my_rates.k10 * gp.H2II_density * gp.HI_density 
		    - 2.f *my_rates.k22 * gp.HI_density * gp.HI_density * gp.HI_density
		    + 0.5 *my_rates.k11 * gp.HII_density * gp.H2I_density 
		    + my_rates.k12 * gp.e_density * gp.H2I_density 
		    + my_rates.k13 * gp.H2I_density * gp.HI_density 
		    + my_rates.k14 * gp.e_density * gp.HM_density 
		    + my_rates.k15 * gp.HM_density * gp.HI_density
		    + 2.f *my_rates.k16 * gp.HII_density * gp.HM_density 
		    + my_rates.k18 * gp.e_density * gp.H2II_density
		    + 0.5 * my_rates.k19 * gp.HM_density * gp.H2II_density;
	    if (chemistry->UVbackground > 0) {
		HIdot += my_rates.k31shield * gp.H2I_density;
	    }
	    if(chemistry->use_dust_evol && gp.dust2gas > tiny) {
		HIdot -= 2.f *my_rates.h2dust * (gp.HI_density + gp.HII_density);
	    }
	}

	return HIdot;
}


void evolve_helium(grackle_part_data *p, grackle_part_data *gp_old, chemistry_data *chemistry, chemistry_rate_storage my_rates, double dtit) /* from step_rate_g() */
{ 
	double scoef, acoef;

	if (!isfinite(p->RT_HeI_ionization_rate)) {
    	fprintf(stderr, "WARNING: Invalid RT_HeI_ionization_rate (%e), setting to 0.\n",
            p->RT_HeI_ionization_rate);
	exit(1);
    	p->RT_HeI_ionization_rate = 0.0; // or a capped value like 1e10
	}

	if (!isfinite(p->RT_HeII_ionization_rate)) {
    	fprintf(stderr, "WARNING: Invalid RT_HeII_ionization_rate (%e), setting to 0.\n",
            p->RT_HeII_ionization_rate);
	exit(1);
    	p->RT_HeII_ionization_rate = 0.0; // or a capped value like 1e10
	}
	
	//printf("my_rates.k26 = %e\n", my_rates.k26);

	/* HeI */
	//scoef = my_rates.k4 * p->HeII_density * p->e_density;
	//acoef = my_rates.k3 * p->e_density + my_rates.k26;
	//if (chemistry->use_radiative_transfer) acoef += p->RT_HeI_ionization_rate;
	double term1 = my_rates.k4 * p->HeII_density;
    	double term2 = term1 * p->e_density;
    	scoef = term2;

    	double term3 = my_rates.k3 * p->e_density;
    	
      if (chemistry->UVbackground > 0) {
          double term4 = term3 + my_rates.k26;
    	  acoef = term4;
	  if (isnan(my_rates.k26) && chemistry->UVbackground > 0) fprintf(stderr, "Error: NaN detected in my_rates.k26 = %g\n", my_rates.k26);
	  if (isnan(term4) && chemistry->UVbackground > 0) fprintf(stderr, "Error: NaN detected in term4 (term3 + k26) = %g\n", term4);
      }

    	if (chemistry->use_radiative_transfer) {
        	double term5 = p->RT_HeI_ionization_rate;
        	acoef += term5;

        	if (isnan(term5)) fprintf(stderr, "Error: NaN detected in RT_HeI_ionization_rate = %g\n", term5);
		//if (isnan(term5)) printf("Thread %d: got NaN\n", omp_get_thread_num());
    	}

    	/* Check for NaN values */
    	if (isnan(my_rates.k4)) fprintf(stderr, "Error: NaN detected in my_rates.k4 = %g\n", my_rates.k4);
    	if (isnan(p->HeII_density)) fprintf(stderr, "Error: NaN detected in p->HeII_density = %g\n", p->HeII_density);
    	if (isnan(p->e_density)) fprintf(stderr, "Error: NaN detected in p->e_density = %g\n", p->e_density);
    	if (isnan(term1)) fprintf(stderr, "Error: NaN detected in term1 (k4 * HeII_density) = %g\n", term1);
    	if (isnan(term2)) fprintf(stderr, "Error: NaN detected in term2 (term1 * e_density) = %g\n", term2);
    	if (isnan(scoef)) fprintf(stderr, "Error: NaN detected in scoef = %g\n", scoef);

    	if (isnan(my_rates.k3)) fprintf(stderr, "Error: NaN detected in my_rates.k3 = %g\n", my_rates.k3);
    	if (isnan(term3)) fprintf(stderr, "Error: NaN detected in term3 (k3 * e_density) = %g\n", term3);
    	if (isnan(acoef)) fprintf(stderr, "Error: NaN detected in acoef = %g\n", acoef);

		
	if (isnan(scoef) || isnan(acoef)) {
        fprintf(stderr, "Error: NaN detected in HeI calculation!\n  scoef = %g, acoef = %g\n", scoef, acoef);
    	}
	
	double HeIp = (scoef * dtit + p->HeI_density) / (1.f + acoef * dtit);
	if (HeIp < 0. ) HeIp = 0.;
	p->delta_HeI = HeIp - p->HeI_density;

	/* HeII */
	scoef = my_rates.k3 * HeIp * p->e_density +
		my_rates.k6 * p->HeIII_density * p->e_density;
	acoef = (my_rates.k4 + my_rates.k5) * p->e_density;
	if (chemistry->UVbackground > 0) {
	    scoef += my_rates.k26 * HeIp;
	    acoef += my_rates.k25;
	}
	if (chemistry->use_radiative_transfer) {
	    scoef += p->RT_HeI_ionization_rate * HeIp;
	    acoef += p->RT_HeII_ionization_rate;
	}

	if (isnan(scoef) || isnan(acoef)) {
        fprintf(stderr, "Error: NaN detected in HeII calculation!\n  scoef = %g, acoef = %g\n", scoef, acoef);
    	}

	double HeIIp = (scoef * dtit + p->HeII_density) / (1.f + acoef * dtit);
	if (HeIIp < 0. ) HeIIp = 0.;
	p->delta_HeII = HeIIp - p->HeII_density;

	/* HeIII */
	scoef = my_rates.k5 * HeIIp * p->e_density;
	acoef = my_rates.k6 * p->e_density;
	if (chemistry->UVbackground > 0) {
	    scoef += my_rates.k25shield * HeIIp;
	}
	if (chemistry->use_radiative_transfer) scoef += p->RT_HeII_ionization_rate * HeIIp;
	
	if (isnan(scoef) || isnan(acoef)) {
        fprintf(stderr, "Error: NaN detected in HeIII calculation!\n  scoef = %g, acoef = %g\n", scoef, acoef);
    	}

	p->delta_HeIII = (scoef * dtit + p->HeIII_density) / (1.f + acoef * dtit) - p->HeIII_density; 
	if (p->delta_HeIII + p->HeIII_density < 0.f ) p->delta_HeIII = -p->HeIII_density;

	/*if (p->delta_HeII != p->delta_HeII || p->delta_HeIII != p->delta_HeIII) {
	    fprintf(stdout, "TROUBLE: dHeII=%g dHeIII=%g scoef=%g acoef=%g HeIIp=%g e=%g k3=%g k4=%g k5=%g k6=%g k25=%g k25shielf=%g k26=%g\n", p->delta_HeII, p->delta_HeIII, scoef, acoef, HeIIp, p->e_density, my_rates.k3, my_rates.k4, my_rates.k5, my_rates.k6, my_rates.k25, my_rates.k25shield, my_rates.k26);
	    fflush(stdout);
	    assert(p->delta_HeII == p->delta_HeII);
	    assert(p->delta_HeIII == p->delta_HeIII);
	}*/
	
	/* Check for NaN values */
    	if (isnan(p->delta_HeII) || isnan(p->delta_HeIII)) {
        fprintf(stderr, "Error: NaN detected in helium evolution!\n");
        fprintf(stderr, "  delta_HeII = %g\n", p->delta_HeII);
        fprintf(stderr, "  delta_HeIII = %g\n", p->delta_HeIII);
        fprintf(stderr, "  Inputs:\n"
                        "    HeI_density = %g\n"
                        "    HeII_density = %g\n"
                        "    HeIII_density = %g\n"
                        "    e_density = %g\n"
                        "    dtit = %g\n"
                        "    RT_HeI_ionization_rate = %g\n"
                        "    RT_HeII_ionization_rate = %g\n",
                        p->HeI_density, p->HeII_density, p->HeIII_density,
                        p->e_density, dtit,
                        p->RT_HeI_ionization_rate, p->RT_HeII_ionization_rate);
        abort();  // Uncomment if you want to stop execution on error
    	}

	return;
}

void check_hydrogen(grackle_part_data *p) 
{
	double H_density = p->HI_density + p->HII_density + p->HM_density;
	double H2_density = p->H2I_density + p->H2II_density;
	double H_frac = (H_density + H2_density) / p->density;
	
	if (isnan(H_frac) || H_frac <= 0.4 || H_frac >= 0.8) {
        	fprintf(stderr, "Assertion failed in check_hydrogen:\n");
        	fprintf(stderr, "  HI_density  = %g\n", p->HI_density);
        	fprintf(stderr, "  HII_density = %g\n", p->HII_density);
        	fprintf(stderr, "  HM_density  = %g\n", p->HM_density);
        	fprintf(stderr, "  H2I_density = %g\n", p->H2I_density);
        	fprintf(stderr, "  H2II_density = %g\n", p->H2II_density);
        	fprintf(stderr, "  Total density = %g\n", p->density);
        	fprintf(stderr, "  H_frac = %g (expected between 0.4 and 0.8)\n", H_frac);
        	if (isnan(H_frac))
            		fprintf(stderr, "  Warning: H_frac is NaN!\n");
		fflush(stderr);  // Make sure output is flushed before crash
		assert(!isnan(H_frac) && H_frac > 0.4 && H_frac < 0.8);
		//assert(0);  // Triggers abort so you still stop execution
    	}
	assert(H_frac > 0.4 && H_frac < 0.8);
	assert(p->HI_density >= 0.f);
	assert(p->HII_density >= 0.f);
}

void evolve_hydrogen(grackle_part_data *p, grackle_part_data *gp_old, chemistry_data *chemistry, chemistry_rate_storage my_rates, double dtit) /* from step_rate_g() */
{ 
	double scoef, acoef;
	
	if (isnan(p->HI_density) || isnan(p->HII_density) || isnan(p->HM_density)) {
    		printf("About to call check_hydrogen in evolve_hydrogen: found NaN in particle!\n");
    		printf("  HI_density  = %e\n", p->HI_density);
		printf("  HII_density = %e\n", p->HII_density);
		printf("  HM_density  = %e\n", p->HM_density);
		fflush(stdout);
	}
	check_hydrogen(p);
	
	if (!isfinite(p->RT_HI_ionization_rate)) {
    	fprintf(stderr, "WARNING: Invalid RT_HI_ionization_rate (%e), setting to 0.\n",
            p->RT_HI_ionization_rate);
	exit(1);
    	p->RT_HI_ionization_rate = 0.0;
	}
	/* HI */
	scoef = my_rates.k2 * p->HII_density * p->e_density;
	if (chemistry->primordial_chemistry >= 2) {
	    scoef += my_rates.k13 * p->HI_density * p->H2I_density
     	          + 0.5 * my_rates.k11 * p->HII_density * p->H2I_density
     	          + my_rates.k12 * p->e_density * p->H2I_density
     	          + my_rates.k14 * p->HM_density * p->e_density
     	          + my_rates.k15 * p->HM_density * p->HI_density
     	          + 2. * my_rates.k16 * p->HM_density * p->HII_density
     	          + my_rates.k18 * p->H2II_density * p->e_density
     	          + 0.5 * my_rates.k19 * p->H2II_density * p->HM_density;
	    if (chemistry->UVbackground > 0) {
     	          scoef += my_rates.k31shield * p->H2I_density;
	    }
	}
	acoef = my_rates.k1 * p->e_density
	      +	my_rates.k57 * p->HI_density
	      +	0.25 * my_rates.k58 * p->HeI_density;

	if (chemistry->UVbackground > 0) {
     	      acoef += my_rates.k24shield;
	}
	if (chemistry->use_radiative_transfer) {
	      acoef += p->RT_HI_ionization_rate;
	}
	if (chemistry->primordial_chemistry >= 2) {
	    acoef += my_rates.k7 * p->e_density
	          + my_rates.k8 * p->HM_density
	          + my_rates.k9 * p->HII_density
	          + 0.5 * my_rates.k10 * p->H2II_density
	          + 2. * my_rates.k22 * p->HI_density * p->HI_density;
	}
	if (chemistry->use_dust_evol  && p->dust2gas > tiny) {
	    acoef += 2. * my_rates.h2dust * p->rhoH;
	}
	double HIp = (scoef * dtit + p->HI_density) / (1.f + acoef * dtit);
	if (HIp < 0. ) HIp = 0.;
	p->delta_HI = HIp - p->HI_density;

	/* HII */
	scoef = my_rates.k1 * p->HI_density * p->e_density
	      +	0.5 * my_rates.k10 * p->H2II_density * p->HI_density
	      +	my_rates.k57 * p->HI_density * p->HI_density
	      +	0.25 * my_rates.k58 * p->HeI_density * p->HI_density;
	acoef = my_rates.k2 * p->e_density
	      +	my_rates.k9 * p->HI_density
	      +	0.5 * my_rates.k11 * p->H2I_density
	      +	my_rates.k16 * p->HM_density
	      +	my_rates.k17 * p->HM_density;

	if (chemistry->UVbackground > 0) {
	    scoef += my_rates.k24shield * p->HI_density;
	}
	if (chemistry->use_radiative_transfer) {
	    scoef += p->RT_HI_ionization_rate * p->HI_density;
	}
	double HIIp = (scoef * dtit + p->HII_density) / (1.f + acoef * dtit);
	if (HIIp < 0. ) HIIp = 0.;
	p->delta_HII = HIIp - p->HII_density;

	if (isnan(p->delta_HII)) {
    		printf("ERROR: p->delta_HII is NaN!\n");

    		// Print densities
    		printf("HI_density = %e, HII_density = %e, e_density = %e\n",
           	p->HI_density, p->HII_density, p->e_density);
    		printf("H2I_density = %e, H2II_density = %e, HM_density = %e, HeI_density = %e\n",
           	p->H2I_density, p->H2II_density, p->HM_density, p->HeI_density);
    		printf("RT_HI_ionization_rate = %e\n", p->RT_HI_ionization_rate);

    		// Print rate coefficients
    		printf("k1 = %e, k2 = %e, k9 = %e, k10 = %e, k11 = %e\n",
           	my_rates.k1, my_rates.k2, my_rates.k9, my_rates.k10, my_rates.k11);
    		printf("k16 = %e, k17 = %e, k24shield = %e, k57 = %e, k58 = %e\n",
           	my_rates.k16, my_rates.k17, my_rates.k24shield, my_rates.k57, my_rates.k58);

    		// Print computed values
    		printf("scoef = %e, acoef = %e, dtit = %e\n", scoef, acoef, dtit);
    		double denominator = 1.0 + acoef * dtit;
    		printf("Denominator (1 + acoef * dtit) = %e\n", denominator);
    		printf("Numerator (scoef * dtit + p->HII_density) = %e\n", scoef * dtit + p->HII_density);
	}

	/* electrons */
	scoef = my_rates.k8 * p->HM_density * p->HI_density
	      + my_rates.k15 * p->HM_density * p->HI_density
	      + my_rates.k17 * p->HM_density * p->HII_density
	      +	my_rates.k10 * p->H2II_density * p->HI_density
	      +	my_rates.k57 * p->HI_density * p->HI_density
	      +	0.25 * my_rates.k58 * p->HeI_density * p->HI_density;

	if (chemistry->UVbackground > 0) {
	    scoef += my_rates.k24shield * HIp
     	          + 0.25 * my_rates.k25shield * (p->HeII_density+p->delta_HeII)
     	          + 0.25 * my_rates.k26shield * (p->HeI_density+p->delta_HeI);
	}
	if (chemistry->use_radiative_transfer) {
	    scoef += p->RT_HI_ionization_rate * HIp
	          + 0.25 * p->RT_HeI_ionization_rate * (p->HeI_density+p->delta_HeI)
	          + 0.25 * p->RT_HeII_ionization_rate * (p->HeII_density+p->delta_HeII);
	}
	acoef = my_rates.k1 * p->HI_density
	      - my_rates.k2 * p->HII_density
	      +	0.25 * my_rates.k3 * p->HeI_density
	      -	0.25 * my_rates.k6 * p->HeIII_density
	      +	0.25 * my_rates.k5 * p->HeII_density
	      -	0.25 * my_rates.k4 * p->HeII_density
	      +	my_rates.k14 * p->HM_density
	      -	my_rates.k7 * p->HI_density
	      -	0.5 * my_rates.k18 * p->H2II_density;
	acoef *= -1.;
	p->delta_e = (scoef * dtit + p->e_density) / (1.f + acoef * dtit) - p->e_density;

	return;
}

void evolve_H2(grackle_part_data *p, int ism_flag, chemistry_data *chemistry, chemistry_rate_storage my_rates, crackle_units cunits, double dtit) /* from step_rate_g() */
{ 
	double scoef, acoef;

	/* H2 immediately destroyed outside ISM and returned to atomic H */
	if (ism_flag == 0) {
	    if (p->rhoH2 > 0.) {
	        p->HI_density += p->H2I_density + p->HM_density;
	        p->HII_density += p->H2II_density;
	        p->H2I_density = p->H2II_density = p->HM_density = 0.;  
	        compute_electron_density(p);
	    }
	    return;
	    } else {
    //printf("ism_flag = %d, rhoH2 = %e, HI_density = %e, HII_density = %e\n",
      //     ism_flag, p->rhoH2, p->HI_density, p->HII_density);
	}

	/* H2 */
	scoef = 2.f * my_rates.k8 * p->HM_density * p->HI_density +
		my_rates.k10 * p->H2II_density * p->HI_density +
		my_rates.k19 * p->H2II_density * p->HM_density +
		2.f * my_rates.k22 * p->HI_density * p->HI_density * p->HI_density;
	acoef = my_rates.k13 * p->HI_density + 
		my_rates.k11 * p->HII_density +
		my_rates.k12 * p->e_density;
	if (chemistry->UVbackground > 0) {
	    acoef += my_rates.k29shield + my_rates.k31shield;
	}
	if (chemistry->use_dust_evol && p->dust2gas > tiny) {
	    scoef += 2.f * my_rates.h2dust * p->HI_density * p->rhoH;
	}
	p->delta_H2I = (scoef * dtit + p->H2I_density) / (1.f + acoef * dtit) - p->H2I_density;
	if (p->delta_H2I < -p->H2I_density) p->delta_H2I = -p->H2I_density;
	if (p->verbose) printf("H2: %g %g %g %g %g %g %g\n",p->density,scoef*dtit,acoef*dtit,my_rates.h2dust,my_rates.k31shield,p->H2I_density,p->delta_H2I);

	/* H- */
	scoef = my_rates.k7 * p->e_density * p->HI_density;
	acoef = (my_rates.k8 + my_rates.k15) * p->HI_density + 
		(my_rates.k16 + my_rates.k17) * p->HII_density + 
		my_rates.k14 * p->e_density +
		0.5 * my_rates.k19 * p->H2II_density;
	if (chemistry->UVbackground > 0) {
	    acoef += my_rates.k27;
	}
	p->delta_HM = (scoef * dtit + p->HM_density) / (1.f + acoef * dtit) - p->HM_density;

	/* H2+ */
	/* For some reason H2+ uses the updated densities, so let's calculate those */
	double HIp, HIIp, H2Ip, HMp, dep;
	HIp = p->HI_density + p->delta_HI;
	HIIp = p->HII_density + p->delta_HII;
	H2Ip = p->H2I_density + p->delta_H2I;
	HMp = p->HM_density + p->delta_HM;
	dep = p->e_density + p->delta_e;
	if (chemistry->UVbackground > 0) {
	    p->delta_H2II = 2.f * (my_rates.k9 * HIp * HIIp +
		0.5 * my_rates.k11 * H2Ip * HIIp +
		my_rates.k17 * HMp * HIIp +
		my_rates.k29shield * H2Ip) / 
		(my_rates.k10 * HIp + my_rates.k18 * dep +
		my_rates.k19 * HMp + my_rates.k28shield + my_rates.k30shield)
		- p->H2II_density;
	}
	else {
	    p->delta_H2II = 2.f * (my_rates.k9 * HIp * HIIp +
		0.5 * my_rates.k11 * H2Ip * HIIp +
		my_rates.k17 * HMp * HIIp) / 
		(my_rates.k10 * HIp + my_rates.k18 * dep +
		my_rates.k19 * HMp)
		- p->H2II_density;
  }
  
	if (isnan(p->delta_H2II)) {
   		printf("\n==== NaN DETECTED in delta_H2II! ====\n");
    		printf("Inputs to delta_H2II calculation:\n");

    		printf("HIp  = %e\n", HIp);
    		printf("HIIp = %e\n", HIIp);
    		printf("H2Ip = %e\n", H2Ip);
    		printf("HMp  = %e\n", HMp);
    		printf("dep  = %e\n", dep);
    		printf("k9   = %e\n", my_rates.k9);
    		printf("k11  = %e\n", my_rates.k11);
    		printf("k17  = %e\n", my_rates.k17);
    		printf("k29shield = %e\n", my_rates.k29shield);
    		printf("k10  = %e\n", my_rates.k10);
    		printf("k18  = %e\n", my_rates.k18);
    		printf("k19  = %e\n", my_rates.k19);
    		printf("k28shield = %e\n", my_rates.k28shield);
    		printf("k30shield = %e\n", my_rates.k30shield);
    		printf("H2II_density (before) = %e\n", p->H2II_density);

    		double numerator = 2.f * (my_rates.k9 * HIp * HIIp +
                              0.5 * my_rates.k11 * H2Ip * HIIp +
                              my_rates.k17 * HMp * HIIp +
                              my_rates.k29shield * H2Ip);
    		double denominator = my_rates.k10 * HIp +
                         my_rates.k18 * dep +
                         my_rates.k19 * HMp +
                         my_rates.k28shield +
                         my_rates.k30shield;

    		printf("Numerator   = %e\n", numerator);
   		printf("Denominator = %e\n", denominator);

    		fflush(stdout);  // make sure it's printed even if crash follows
	}

	/*if (p->delta_H2I != p->delta_H2I || p->delta_H2II != p->delta_H2II) {
	    fprintf(stdout, "TROUBLE: dH2I=%g dH2II=%g scoef=%g acoef=%g H2Ip=%g e=%g k3=%g k4=%g k5=%g k6=%g k25=%g k25shielf=%g k26=%g\n", p->delta_H2I, p->delta_H2II, scoef, acoef, H2Ip, p->e_density, my_rates.k3, my_rates.k4, my_rates.k5, my_rates.k6, my_rates.k25, my_rates.k25shield, my_rates.k26);
	    fflush(stdout);
	    assert(p->delta_H2I == p->delta_H2I);
	    assert(p->delta_H2II == p->delta_H2II);
	}*/
	return;
}

void print_nan_check(const grackle_part_data *gp, const char *label) {
    printf("\n=== NaN Check: %s ===\n", label);
    printf("HI_density   = %e\n", gp->HI_density);
    printf("HII_density  = %e\n", gp->HII_density);
    printf("H2I_density  = %e\n", gp->H2I_density);
    printf("H2II_density = %e\n", gp->H2II_density);
    printf("HM_density   = %e\n", gp->HM_density);
    printf("HeI_density  = %e\n", gp->HeI_density);
    printf("HeII_density = %e\n", gp->HeII_density);
    printf("HeIII_density= %e\n", gp->HeIII_density);
    printf("rhoH         = %e\n", gp->rhoH);
    printf("delta_HI     = %e\n", gp->delta_HI);
    printf("delta_HII    = %e\n", gp->delta_HII);
    printf("delta_H2I    = %e\n", gp->delta_H2I);
    printf("delta_H2II   = %e\n", gp->delta_H2II);
    printf("delta_HM     = %e\n", gp->delta_HM);
    printf("delta_HeI    = %e\n", gp->delta_HeI);
    printf("delta_HeII   = %e\n", gp->delta_HeII);
    printf("delta_HeIII  = %e\n", gp->delta_HeIII);
    fflush(stdout);
}

void check_nan_and_print(const grackle_part_data *gp, const char *label) {
    if (
        isnan(gp->HI_density)   || isnan(gp->HII_density)  || isnan(gp->H2I_density) ||
        isnan(gp->H2II_density) || isnan(gp->HM_density)   ||
        isnan(gp->HeI_density)  || isnan(gp->HeII_density) || isnan(gp->HeIII_density) ||
        isnan(gp->rhoH)         ||
        isnan(gp->delta_HI)     || isnan(gp->delta_HII)    || isnan(gp->delta_H2I) ||
        isnan(gp->delta_H2II)   || isnan(gp->delta_HM)     ||
        isnan(gp->delta_HeI)    || isnan(gp->delta_HeII)   || isnan(gp->delta_HeIII)
    ) {
        print_nan_check(gp, label);
    }
}

void evolve_elements(grackle_part_data *gp, grackle_part_data *gp_old, chemistry_data *chemistry)
{
        check_nan_and_print(gp, "Before limiting changes");
	check_hydrogen(gp);
	
	/* Limit change for H species */
	if (gp->HI_density + gp->delta_HI < 0.) gp->delta_HI = -gp->HI_density;
	if (gp->HI_density + gp->delta_HI > gp->rhoH) gp->delta_HI = gp->rhoH - gp->HI_density;
	if (gp->HII_density + gp->delta_HII < 0.) gp->delta_HII = -gp->HII_density;
	if (gp->HII_density + gp->delta_HII > gp->rhoH) gp->delta_HII = gp->rhoH - gp->HII_density;
	if (gp->H2I_density + gp->delta_H2I < 0.) gp->delta_H2I = -gp->H2I_density;
	if (gp->H2I_density + gp->delta_H2I > gp->rhoH) gp->delta_H2I = gp->rhoH - gp->H2I_density;
	if (gp->H2II_density + gp->delta_H2II < 0.) gp->delta_H2II = -gp->H2II_density;
	if (gp->H2II_density + gp->delta_H2II > gp->rhoH) gp->delta_H2II = gp->rhoH - gp->H2II_density;
	if (gp->HM_density + gp->delta_HM < 0.) gp->delta_HM = -gp->HM_density;
	if (gp->HM_density + gp->delta_HM > gp->rhoH) gp->delta_HM = gp->rhoH - gp->HM_density;

	/* Evolve H species */
	gp->HI_density += gp->delta_HI;
	gp->HII_density += gp->delta_HII;
	gp->H2I_density += gp->delta_H2I;
	gp->H2II_density += gp->delta_H2II;
	gp->HM_density += gp->delta_HM;

	/* Evolve HeI, HeII */
	gp->HeI_density += gp->delta_HeI;
	gp->HeII_density += gp->delta_HeII;
	gp->HeIII_density += gp->delta_HeIII;
	if (gp->HeI_density < tiny) gp->HeI_density = tiny;
	if (gp->HeII_density < tiny) gp->HeII_density = tiny;
	if (gp->HeIII_density < tiny) gp->HeIII_density = tiny;

	/* Ensure H & He are conserved (since these are unaffected by dust) */
	const double H_old = gp_old->HI_density + gp_old->HII_density + gp_old->H2I_density + gp_old->H2II_density + gp_old->HM_density;
	const double He_old = gp_old->HeI_density + gp_old->HeII_density + gp_old->HeIII_density;
	const double H_new = gp->HI_density + gp->HII_density + gp->H2I_density + gp->H2II_density + gp->HM_density;
	const double He_new = gp->HeI_density + gp->HeII_density + gp->HeIII_density;

	const double H_factor = H_old / (H_new+tiny);
	const double He_factor = He_old / (He_new+tiny);

	/* Correct densities to conserve H, He */
	gp->HI_density *= H_factor;
	gp->HII_density *= H_factor;
	gp->H2I_density *= H_factor;
	gp->H2II_density *= H_factor;
	gp->HM_density *= H_factor;

	gp->HeI_density *= He_factor;
	gp->HeII_density *= He_factor;
	gp->HeIII_density *= He_factor;
	
	check_nan_and_print(gp, "After conservation correction");
	check_hydrogen(gp);

	/* Compute new electron density */
	compute_electron_density(gp);

	return;
}

void evolve_pred_corr(grackle_part_data *gp, grackle_part_data *gp_old, chemistry_data *chemistry)
{
	/* Evolve all species */
	gp->HI_density = 0.5 * (gp->HI_density + gp_old->HI_density);
	gp->HII_density = 0.5 * (gp->HII_density + gp_old->HII_density);
	gp->HeI_density = 0.5 * (gp->HeI_density + gp_old->HeI_density);
	gp->HeII_density = 0.5 * (gp->HeII_density + gp_old->HeII_density);
	gp->HeIII_density = 0.5 * (gp->HeIII_density + gp_old->HeIII_density);

	/* Now ensure H & He are conserved */
	const double H_old = gp_old->HI_density + gp_old->HII_density + gp_old->H2I_density + gp_old->H2II_density + gp_old->HM_density;
	const double He_old = gp_old->HeI_density + gp_old->HeII_density + gp_old->HeIII_density;
	const double H_new = gp->HI_density + gp->HII_density + gp->H2I_density + gp->H2II_density + gp->HM_density;
	const double He_new = gp->HeI_density + gp->HeII_density + gp->HeIII_density;

	const double H_factor = H_old / (H_new+tiny);
	const double He_factor = He_old / (He_new+tiny);

	/* Correct densities proportionally to conserve H, He */
	gp->HI_density *= H_factor;
	gp->HII_density *= H_factor;
	gp->H2I_density *= H_factor;
	gp->H2II_density *= H_factor;
	gp->HM_density *= H_factor;

	gp->HeI_density *= He_factor;
	gp->HeII_density *= He_factor;
	gp->HeIII_density *= He_factor;

	/* Compute new electron density */
	compute_electron_density(gp);

	return;
}

void crackle_cooling_time(grackle_field_data *p, chemistry_data *chemistry, chemistry_data_storage grackle_rates, photo_rate_storage my_uvb_rates, code_units *units, gr_float *tcool)
{
	grackle_part_data gp[1];
	chemistry_rate_storage my_rates;  // interpolated rates for this field
	interp_struct interpolation;
	crackle_units cunits;

	/* Copy to grackle_part_data */
	copy_grackle_fields_to_part(p, gp, chemistry);
	/* Compute some basic properties */
	set_rhot(gp, units, chemistry);
	/* Set up various unit conversions etc */
	set_crackle_units(units, grackle_rates, chemistry->Gamma, &cunits);
	/* initialize inteprolation of chemistry rate tables for this particle */
	init_temperature_interpolation(gp, chemistry, &interpolation, cunits, grackle_rates);
	/* Set up cooling/heating rates interpolation */
	setup_temperature_interpolation(gp->tgas, chemistry, &interpolation);
	/* Get interpolated chemistry rates for this particle */
	lookup_chemistry_coeffs(chemistry->primordial_chemistry, grackle_rates, &my_rates, &interpolation);  
	/* ISM flag: see whether we should be forming dust and evolving H2 */
	int ism_flag = (gp->isrf_habing >= 0.);
	if (gp->isrf_habing < 0.) gp->isrf_habing = 0.;
	/* Compute rate of change of thermal energy */
	gp->verbose = 0;
	compute_edot(gp, chemistry, grackle_rates, &my_rates, my_uvb_rates, &interpolation, units, cunits);  
	gp->verbose = 0;

	*tcool = gp->internal_energy * gp->density / (gp->edot+tiny);
	return;
}

void crackle_temperature(grackle_field_data *p, chemistry_data *chemistry, code_units *units, gr_float *temp)
{
	grackle_part_data gp;

	/* Copy to grackle_part_data */
	copy_grackle_fields_to_part(p, &gp, chemistry);
	/* Compute some basic properties */
	set_rhot(&gp, units, chemistry);

	*temp = gp.tgas;
	return;
}

