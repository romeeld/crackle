

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>

bool has_nan_in_gp(const grackle_part_data *gp) {
    bool found_nan = false;

    // Macro to reduce repetition
    #define CHECK_NAN(field) \
        if (isnan(gp->field)) { \
            printf("NaN detected in gp->%s\n", #field); \
            found_nan = true; \
        }

    CHECK_NAN(density);
    CHECK_NAN(HI_density);
    CHECK_NAN(HII_density);
    CHECK_NAN(HM_density);
    CHECK_NAN(HeI_density);
    CHECK_NAN(HeII_density);
    CHECK_NAN(HeIII_density);
    CHECK_NAN(H2I_density);
    CHECK_NAN(H2II_density);
    CHECK_NAN(DI_density);
    CHECK_NAN(DII_density);
    CHECK_NAN(HDI_density);
    CHECK_NAN(e_density);
    CHECK_NAN(metal_density);
    CHECK_NAN(internal_energy);
    CHECK_NAN(x_velocity);
    CHECK_NAN(y_velocity);
    CHECK_NAN(z_velocity);
    CHECK_NAN(volumetric_heating_rate);
    CHECK_NAN(specific_heating_rate);
    CHECK_NAN(temperature_floor);
    CHECK_NAN(RT_heating_rate);
    CHECK_NAN(RT_HI_ionization_rate);
    CHECK_NAN(RT_HeI_ionization_rate);
    CHECK_NAN(RT_HeII_ionization_rate);
    CHECK_NAN(RT_H2_dissociation_rate);
    CHECK_NAN(H2_self_shielding_length);
    CHECK_NAN(H2_custom_shielding_factor);
    CHECK_NAN(isrf_habing);
    CHECK_NAN(dust_density);
    CHECK_NAN(SNe_density);
    CHECK_NAN(tdust);

    CHECK_NAN(edot);
    CHECK_NAN(edot_ext);
    CHECK_NAN(dedot);
    CHECK_NAN(HIdot);
    CHECK_NAN(fSShHI);
    CHECK_NAN(fSShHeI);
    CHECK_NAN(fSShHeII);
    CHECK_NAN(rhoH);
    CHECK_NAN(rhoHe);
    CHECK_NAN(rhoH2);
    CHECK_NAN(mmw);
    CHECK_NAN(tgas);
    CHECK_NAN(logtem);
    CHECK_NAN(delta_HI);
    CHECK_NAN(delta_HII);
    CHECK_NAN(delta_HeI);
    CHECK_NAN(delta_HeII);
    CHECK_NAN(delta_HeIII);
    CHECK_NAN(delta_HM);
    CHECK_NAN(delta_H2I);
    CHECK_NAN(delta_H2II);
    CHECK_NAN(delta_e);

    for (int i = 0; i < NUM_METAL_SPECIES_GRACKLE; i++) {
        if (isnan(gp->gas_metalDensity[i])) {
            printf("NaN detected in gp->gas_metalDensity[%d]\n", i);
            found_nan = true;
        }
        if (isnan(gp->dust_metalDensity[i])) {
            printf("NaN detected in gp->dust_metalDensity[%d]\n", i);
            found_nan = true;
        }
    }

    return found_nan;

    #undef CHECK_NAN
}

int fix_nan_in_gp(grackle_part_data *gp, FILE *logfile) {
    int count = 0;

    #define FIX_NAN(field) \
        if (isnan(gp->field)) { \
            gp->field = 0.0; \
            count++; \
            if (logfile) fprintf(logfile, "NaN fixed in gp->%s\n", #field); \
        }

    FIX_NAN(density);
    FIX_NAN(HI_density);
    FIX_NAN(HII_density);
    FIX_NAN(HM_density);
    FIX_NAN(HeI_density);
    FIX_NAN(HeII_density);
    FIX_NAN(HeIII_density);
    FIX_NAN(H2I_density);
    FIX_NAN(H2II_density);
    FIX_NAN(DI_density);
    FIX_NAN(DII_density);
    FIX_NAN(HDI_density);
    FIX_NAN(e_density);
    FIX_NAN(metal_density);
    FIX_NAN(internal_energy);
    FIX_NAN(x_velocity);
    FIX_NAN(y_velocity);
    FIX_NAN(z_velocity);
    FIX_NAN(volumetric_heating_rate);
    FIX_NAN(specific_heating_rate);
    FIX_NAN(temperature_floor);
    FIX_NAN(RT_heating_rate);
    FIX_NAN(RT_HI_ionization_rate);
    FIX_NAN(RT_HeI_ionization_rate);
    FIX_NAN(RT_HeII_ionization_rate);
    FIX_NAN(RT_H2_dissociation_rate);
    FIX_NAN(H2_self_shielding_length);
    FIX_NAN(H2_custom_shielding_factor);
    FIX_NAN(isrf_habing);
    FIX_NAN(dust_density);
    FIX_NAN(SNe_density);
    FIX_NAN(tdust);

    FIX_NAN(edot);
    FIX_NAN(edot_ext);
    FIX_NAN(dedot);
    FIX_NAN(HIdot);
    FIX_NAN(fSShHI);
    FIX_NAN(fSShHeI);
    FIX_NAN(fSShHeII);
    FIX_NAN(rhoH);
    FIX_NAN(rhoHe);
    FIX_NAN(rhoH2);
    FIX_NAN(mmw);
    FIX_NAN(tgas);
    FIX_NAN(logtem);
    FIX_NAN(delta_HI);
    FIX_NAN(delta_HII);
    FIX_NAN(delta_HeI);
    FIX_NAN(delta_HeII);
    FIX_NAN(delta_HeIII);
    FIX_NAN(delta_HM);
    FIX_NAN(delta_H2I);
    FIX_NAN(delta_H2II);
    FIX_NAN(delta_e);

    for (int i = 0; i < NUM_METAL_SPECIES_GRACKLE; i++) {
        if (isnan(gp->gas_metalDensity[i])) {
            gp->gas_metalDensity[i] = 0.0;
            count++;
            if (logfile) fprintf(logfile, "NaN fixed in gp->gas_metalDensity[%d]\n", i);
        }
        if (isnan(gp->dust_metalDensity[i])) {
            gp->dust_metalDensity[i] = 0.0;
            count++;
            if (logfile) fprintf(logfile, "NaN fixed in gp->dust_metalDensity[%d]\n", i);
        }
    }

    if (count > 0 && logfile) {
        fprintf(logfile, "Total NaNs fixed in this gp structure: %d\n", count);
    }

    return count;
}

static inline void set_crackle_units(code_units *units, chemistry_data_storage grackle_rates, double gamma, crackle_units *cunits) {
                /* Set units */
        cunits->dom      = units->density_units/mh; // converts to equivalent physical H number density
        cunits->dom_inv  = 1./cunits->dom;
        cunits->time_to_cgs = units->time_units;  // tbase1 in original grackle
        cunits->length_to_cgs = units->length_units;  // xbase1
        cunits->density_to_cgs = units->density_units;  // dbase1
        cunits->temp_to_K = units->temperature_units;  // usually 1
        cunits->coolunit = (pow(units->a_units, 5) * units->length_units*units->length_units * mh*mh) / (units->time_units*units->time_units*units->time_units * units->density_units);
        cunits->coolunit_inv = 1./cunits->coolunit;
        cunits->velunit = units->length_units / (units->a_value * units->time_units);
        cunits->chunit = 1.60218e-12 / (2 * cunits->velunit * cunits->velunit * mh);

        /* Compton cooling coeffs */
        cunits->redshift = 1./units->a_value - 1.;
        cunits->zr1      = cunits->redshift + 1.;
        cunits->compton1 = grackle_rates.comp * cunits->zr1*cunits->zr1*cunits->zr1*cunits->zr1;
        cunits->compton2 = TCMB0 * cunits->zr1;
        cunits->log10_TCMB = log10(cunits->compton2);
        cunits->c_ljeans = sqrt((gamma * M_PI * kboltz) / (GravConst * mh * cunits->density_to_cgs));

        return;
}

static inline void compute_electron_density(grackle_part_data *gp)
{
        gp->e_density = gp->HII_density + 0.25 * gp->HeII_density + 0.5 * gp->HeIII_density;
        gp->e_density += 0.5 * gp->H2II_density - gp->HM_density;
	if (gp->e_density < 0.) {  // lower HM_density to make it zero, add to HI
	    gp->HI_density -= gp->e_density;
	    gp->HM_density += gp->e_density;
	    gp->e_density = 0.;
	}
	if (gp->e_density < 0.) {
        fprintf(stderr, "Assertion failed: gp->e_density = %g\n", gp->e_density);
    }
	if (isnan(gp->e_density)) {
    	fprintf(stderr, "Error: e_density is NaN! Check inputs.\n");
	fprintf(stderr, "Assertion failed: gp->HII_density = %g\n"
			"gp->HeII_density = %g\n"
			"gp->HeIII_density = %g\n"
			"gp->H2II_density = %g\n"
			"gp->HM_density = %g\n", gp->HII_density, gp->HeII_density, gp->HeIII_density, gp->H2II_density, gp->HM_density);
    	//abort();  // Or handle the issue in another way
}
	assert(gp->e_density >= 0.);
}

static inline void set_rhot(grackle_part_data *gp, code_units *units, chemistry_data *chemistry)
{
	/* Compute densities and store */
	gp->rhoHe = gp->HeI_density + gp->HeII_density + gp->HeIII_density;
	gp->rhoH2 = gp->H2I_density + gp->H2II_density;
	gp->rhoH = gp->HI_density + gp->HII_density + gp->rhoH2 + gp->HM_density;
	gp->nH = gp->rhoH * units->density_units / mh;

	/* Compute mean molecular weighted density */
	gp->mmw = gp->HI_density + gp->HII_density + 0.25 * gp->rhoHe + gp->e_density;
	gp->mmw += 0.5 * gp->rhoH2 + gp->HM_density;
	gp->mmw += gp->metal_density / 16.;  
	gp->mmw = gp->density / gp->mmw;

	/* Calculate solar-scaled total metallicity */
	gp->metallicity = gp->metal_density / (gp->density * chemistry->SolarMetalFractionByMass);
	if (chemistry->use_dust_evol) gp->dust2gas = gp->dust_density / gp->density;
	else gp->dust2gas = 0.009387 * gp->metallicity; // value from local mol clouds, for metallicity in solar units

	/* Compute temperature */
	gp->tgas = fmin(fmax((chemistry->Gamma - 1.f) * gp->internal_energy * gp->mmw * units->temperature_units, chemistry->TemperatureStart), chemistry->TemperatureEnd);

	/* Correct temperature for H2 */
	const double Gamma1 = chemistry->Gamma - 1.;
	double Gamma2 = 2.5;
	double nH2 = 0.5 * gp->rhoH2;
	double nother = 0.25 * gp->rhoHe + gp->HI_density + gp->HII_density + gp->e_density;
	if (nH2 > 1.e-3 * nother){
	    double gfac = fmin(6100. / gp->tgas, 10.);
	    if (gfac < 10.) {
	        Gamma2 = 0.5 * (5.+ 2.*gfac*gfac*exp(gfac)/((exp(gfac)-1.)*(exp(gfac)-1.)));
	    }
	}
	Gamma2 = (nH2+nother) / (nH2*Gamma2 + nother/Gamma1);
	gp->tgas *= Gamma2 / Gamma1; 
	gp->logtem = log(gp->tgas);
	//printf("rhot: rhoH=%g u=%g mmw=%g tunit=%g tgas=%g t=%g\n",gp->rhoH, gp->internal_energy, gp->mmw, units->temperature_units, (chemistry->Gamma - 1.f) * gp->internal_energy * gp->mmw * units->temperature_units, gp->tgas);
	assert(gp->tgas > 0.);
	compute_electron_density(gp);
}

static inline int apply_temperature_bounds(grackle_part_data *gp, chemistry_data *chemistry, float tfloor, float tceiling)
{
	/* If we've reached temp floor and are still cooling, then return 1 to signal that iteration should end */
        if (gp->tgas < tfloor && gp->edot < 0.) {
            gp->internal_energy *= tfloor / gp->tgas;
            gp->tgas = tfloor; 
            return 1;
        }
        /* If we've reached temp ceiling and are still heating, we're also done */
	else if (gp->tgas > tceiling && gp->edot > 0.) {
            gp->internal_energy *= tceiling / gp->tgas;
            gp->tgas = tceiling;
            return 1;
        }
	return 0;
}

static inline void copy_grackle_fields_to_part(grackle_field_data *p, grackle_part_data *gp, chemistry_data *chemistry) 
{

	gp->grid_rank = p->grid_rank;
	gp->grid_dimension = p->grid_dimension[0];
	gp->grid_start = p->grid_start[0];
	gp->grid_end = p->grid_end[0];

	gp->density = p->density[0];
	gp->HI_density = p->HI_density[0];
	gp->HII_density = p->HII_density[0];
	gp->HM_density = p->HM_density[0];
	gp->HeI_density = p->HeI_density[0];
	gp->HeII_density = p->HeII_density[0];
	gp->HeIII_density = p->HeIII_density[0];
	if (chemistry->primordial_chemistry >= 2) {
	    gp->H2I_density = p->H2I_density[0];
	    gp->H2II_density = p->H2II_density[0];
	}
	else {
	    gp->H2I_density = 0.;
	    gp->H2II_density = 0.;
	}
	if (chemistry->primordial_chemistry >= 3) {
	    printf("Crackle does not currently work with primordial_chemistry >=3, sorry!\n");
	    return;
	    //gp->DI_density = p->DI_density[0];
	    //gp->DII_density = p->DII_density[0];
	    //gp->HDI_density = p->HDI_density[0];
	}
	else {
	    gp->DI_density = 0.;
	    gp->DII_density = 0.;
	    gp->HDI_density = 0.;
	}
	gp->e_density = p->e_density[0];
	gp->metal_density = p->metal_density[0];
	gp->internal_energy = p->internal_energy[0];
	gp->x_velocity = p->x_velocity[0];
	gp->y_velocity = p->y_velocity[0];
	gp->z_velocity = p->z_velocity[0];
	if (chemistry->use_volumetric_heating_rate) gp->volumetric_heating_rate = p->volumetric_heating_rate[0];
	else gp->volumetric_heating_rate = 0.;
	if (chemistry->use_specific_heating_rate) gp->specific_heating_rate = p->specific_heating_rate[0];
	else gp->specific_heating_rate = 0.;
	if (p->temperature_floor != NULL) gp->temperature_floor = p->temperature_floor[0];
	else gp->temperature_floor = tiny;
	if (chemistry->use_radiative_transfer) {
	    gp->RT_heating_rate = p->RT_heating_rate[0];
	    gp->RT_HI_ionization_rate = p->RT_HI_ionization_rate[0];
	    gp->RT_HeI_ionization_rate = p->RT_HeI_ionization_rate[0];
	    gp->RT_HeII_ionization_rate = p->RT_HeII_ionization_rate[0];
	    gp->RT_H2_dissociation_rate = p->RT_H2_dissociation_rate[0];
	}
	else {
	    gp->RT_heating_rate = 0.;
	    gp->RT_HI_ionization_rate = 0.;
	    gp->RT_HeI_ionization_rate = 0.;
	    gp->RT_HeII_ionization_rate = 0.;
	    gp->RT_H2_dissociation_rate = 0.;
	}
	if (chemistry->H2_custom_shielding==2) {
	    gp->H2_self_shielding_length = p->H2_self_shielding_length[0];
	}
	else if (chemistry->H2_custom_shielding==1) {
	    gp->H2_custom_shielding_factor = p->H2_custom_shielding_factor[0];
	}
	else {
	    gp->H2_self_shielding_length = 0.;
	    gp->H2_custom_shielding_factor = 0.;
	}
	if (chemistry->use_isrf_field) gp->isrf_habing = p->isrf_habing[0];
	else gp->isrf_habing = 0.;

	/* Dust model */
	if (chemistry->use_dust_evol) {
	    gp->dust_density = p->dust_density[0];
	    gp->gas_metalDensity[0] = p->He_gas_metalDensity[0];
	    gp->gas_metalDensity[1] = p->C_gas_metalDensity[0];
	    gp->gas_metalDensity[2] = p->N_gas_metalDensity[0];
	    gp->gas_metalDensity[3] = p->O_gas_metalDensity[0];
	    gp->gas_metalDensity[4] = p->Ne_gas_metalDensity[0];
	    gp->gas_metalDensity[5] = p->Mg_gas_metalDensity[0];
	    gp->gas_metalDensity[6] = p->Si_gas_metalDensity[0];
	    gp->gas_metalDensity[7] = p->S_gas_metalDensity[0];
	    gp->gas_metalDensity[8] = p->Ca_gas_metalDensity[0];
	    gp->gas_metalDensity[9] = p->Fe_gas_metalDensity[0];
	    gp->dust_metalDensity[0] = p->He_dust_metalDensity[0];
	    gp->dust_metalDensity[1] = p->C_dust_metalDensity[0];
	    gp->dust_metalDensity[2] = p->N_dust_metalDensity[0];
	    gp->dust_metalDensity[3] = p->O_dust_metalDensity[0];
	    gp->dust_metalDensity[4] = p->Ne_dust_metalDensity[0];
	    gp->dust_metalDensity[5] = p->Mg_dust_metalDensity[0];
	    gp->dust_metalDensity[6] = p->Si_dust_metalDensity[0];
	    gp->dust_metalDensity[7] = p->S_dust_metalDensity[0];
	    gp->dust_metalDensity[8] = p->Ca_dust_metalDensity[0];
	    gp->dust_metalDensity[9] = p->Fe_dust_metalDensity[0];
	    gp->SNe_density = p->SNe_ThisTimeStep[0];
	    assert(gp->dust_density == gp->dust_density);
	}
	else {
	    gp->dust_density = 0.;
	    gp->SNe_density = 0.;
	    for (int i=0; i<NUM_METAL_SPECIES_GRACKLE; i++) {
	    	gp->gas_metalDensity[i] = 0.;
	    	gp->dust_metalDensity[i] = 0.;
	    }
	}

	/* initialize evolved quantities */
	gp->edot = gp->edot_ext = gp->dedot = gp->HIdot = gp->fSShHI = gp->fSShHeI = gp->fSShHeII = gp->rhoH = gp->rhoHe = gp->rhoH2 = gp->mmw = gp->tgas = gp->logtem = gp->delta_HI = gp->delta_HII = gp->delta_HeI = gp->delta_HeII = gp->delta_HeIII = gp->delta_HM = gp->delta_H2I = gp->delta_H2II = gp->delta_e = 0.;
	gp->tdust = -1.;  // initialize to negative value so it makes an initial guess
	gp->verbose = 0;

}

static inline void copy_grackle_fields_from_part(grackle_field_data *p, grackle_part_data *gp, chemistry_data *chemistry) 
{

	p->density[0] = gp->density;
	p->HI_density[0] = gp->HI_density;
	p->HII_density[0] = gp->HII_density;
	p->HM_density[0] = gp->HM_density;
	p->HeI_density[0] = gp->HeI_density;
	p->HeII_density[0] = gp->HeII_density;
	p->HeIII_density[0] = gp->HeIII_density;
	if (chemistry->primordial_chemistry >= 2) {
	    p->H2I_density[0] = gp->H2I_density;;
	    p->H2II_density[0] = gp->H2II_density;
	}
	if (chemistry->primordial_chemistry >= 3) {
	    printf("Crackle does not currently work with primordial_chemistry >=3, sorry!\n");
	    return;
	    //p->DI_density[0] = gp->DI_density;
	    //p->DII_density[0] = gp->DII_density;
	    //p->HDI_density[0] = gp->HDI_density;
	}
	p->e_density[0] = gp->e_density;
	p->metal_density[0] = gp->metal_density;
	p->internal_energy[0] = gp->internal_energy;

	if (chemistry->use_dust_evol) {
	    p->dust_density[0] = gp->dust_density;
	    p->He_gas_metalDensity[0] = gp->gas_metalDensity[0];
	    p->C_gas_metalDensity[0] = gp->gas_metalDensity[1];
	    p->N_gas_metalDensity[0] = gp->gas_metalDensity[2];
	    p->O_gas_metalDensity[0] = gp->gas_metalDensity[3];
	    p->Ne_gas_metalDensity[0] = gp->gas_metalDensity[4];
	    p->Mg_gas_metalDensity[0] = gp->gas_metalDensity[5];
	    p->Si_gas_metalDensity[0] = gp->gas_metalDensity[6];
	    p->S_gas_metalDensity[0] = gp->gas_metalDensity[7];
	    p->Ca_gas_metalDensity[0] = gp->gas_metalDensity[8];
	    p->Fe_gas_metalDensity[0] = gp->gas_metalDensity[9];
	    p->He_dust_metalDensity[0] = gp->dust_metalDensity[0];
	    p->C_dust_metalDensity[0] = gp->dust_metalDensity[1];
	    p->N_dust_metalDensity[0] = gp->dust_metalDensity[2];
	    p->O_dust_metalDensity[0] = gp->dust_metalDensity[3];
	    p->Ne_dust_metalDensity[0] = gp->dust_metalDensity[4];
	    p->Mg_dust_metalDensity[0] = gp->dust_metalDensity[5];
	    p->Si_dust_metalDensity[0] = gp->dust_metalDensity[6];
	    p->S_dust_metalDensity[0] = gp->dust_metalDensity[7];
	    p->Ca_dust_metalDensity[0] = gp->dust_metalDensity[8];
	    p->Fe_dust_metalDensity[0] = gp->dust_metalDensity[9];
	}

}


