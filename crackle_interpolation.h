
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "phys_constants.h"


static inline void compute_self_shielded_rates(grackle_part_data *gp, chemistry_data *chemistry, chemistry_rate_storage *my_rates, photo_rate_storage my_uvb_rates, crackle_units cunits) {

	my_rates->k24 = my_rates->k25 = my_rates->k26 =
		my_rates->k27 = my_rates->k28 = my_rates->k29 =
		my_rates->k30 = my_rates->k31 = my_rates->piHI =
		my_rates->piHeI = my_rates->piHeII = my_rates->crsHI =
		my_rates->crsHeI = my_rates->crsHeII =
		my_rates->comp_xray = my_rates->temp_xray = 0.;

	/* Store everything in my_rates for convenience */
	if (chemistry->UVbackground > 0) {
	    my_rates->k24 = my_uvb_rates.k24;
	    my_rates->k25 = my_uvb_rates.k25;
	    my_rates->k26 = my_uvb_rates.k26;
	    my_rates->k27 = my_uvb_rates.k27;
	    my_rates->k28 = my_uvb_rates.k28;
	    my_rates->k29 = my_uvb_rates.k29;
	    my_rates->k30 = my_uvb_rates.k30;
	    my_rates->k31 = my_uvb_rates.k31;
	    my_rates->piHI = my_uvb_rates.piHI;
	    my_rates->piHeI = my_uvb_rates.piHeI;
	    my_rates->piHeII = my_uvb_rates.piHeII;
	    my_rates->crsHI = my_uvb_rates.crsHI;
	    my_rates->crsHeI = my_uvb_rates.crsHeI;
	    my_rates->crsHeII = my_uvb_rates.crsHeII;
	    my_rates->comp_xray = my_uvb_rates.comp_xray;
	    my_rates->temp_xray = my_uvb_rates.temp_xray;
	}

	/* Set up case for (default) no self-shielding */
        gp->fSShHI = gp->fSShHeI = gp->fSShHeII = 1.f;   
	my_rates->k24shield = my_rates->k24;
	my_rates->k25shield = my_rates->k25;
	my_rates->k26shield = my_rates->k26;
	my_rates->k27shield = my_rates->k27;
	my_rates->k28shield = my_rates->k28;
	my_rates->k29shield = my_rates->k29;
	my_rates->k30shield = my_rates->k30;
	my_rates->k31shield = my_rates->k31;

	/* No self-shielding, so no changes to rates */
	if (chemistry->self_shielding_method == 0 || chemistry->UVbackground == 0) return;

        /* We have self-shielding! */
        /* HI self-shielding factor */
	if (my_rates->k24 < tiny) gp->fSShHI = 1.;
	else {
            const double nSSh =  6.73e-3 * pow(my_uvb_rates.crsHI * cunits.time_to_cgs * 1.e-12 / (2.49e-18 * my_rates->k24), -0.6666667) * pow(gp->tgas*1.e-4, 0.17);
            const double nratio = gp->rhoH * cunits.dom / nSSh;
            gp->fSShHI = 0.98*pow(1.+pow(nratio,1.64), -2.28) + 0.02*pow(1.+nratio, -0.84);
	}

        if (chemistry->self_shielding_method == 2 || chemistry->self_shielding_method == 3) {
            /* HeI self-shielding factor */
	    if (my_rates->k26 < tiny) gp->fSShHeI = 1.;
	    else {
                const double nSSh_he =  6.73e-3 * pow(my_uvb_rates.crsHI * cunits.time_to_cgs * 1.e-12 / (2.49e-18 * my_rates->k26), -0.6666667) * pow(gp->tgas*1.e-4, 0.17);
                const double nratio_he = gp->rhoHe * cunits.dom / nSSh_he;
                gp->fSShHeI = 0.98*pow(1.+pow(nratio_he,1.64), -2.28) + 0.02*pow(1.+nratio_he, -0.84);
	    }
        }

        if (chemistry->self_shielding_method == 3) {
            /* HeII self-shielding factor: in this mode, HeII assumed to be completely shielded */
            gp->fSShHeII = 0.f;
        }

	/* Rahmati+2013 H self-shielding */
	my_rates->k24shield *= gp->fSShHI;
	my_rates->k29shield *= gp->fSShHI;
	/* Rahmati plus assuming He closely follows H */
	if (chemistry->self_shielding_method == 2 || chemistry->self_shielding_method == 3) {
	    my_rates->k26shield *= gp->fSShHeI;
	    my_rates->k28shield *= gp->fSShHeI;
	    my_rates->k30shield *= gp->fSShHeI;
	}

	/* Set up H2 self-shielding */
	double l_H2shield;
	if (chemistry->H2_self_shielding > 0 && chemistry->H2_custom_shielding != 1) {
	    if (chemistry->H2_custom_shielding > 0) {
	        l_H2shield = gp->H2_self_shielding_length * cunits.length_to_cgs; // user specifies the H2 shielding length
	    }
	    else {  
	        l_H2shield = cunits.c_ljeans * sqrt(gp->tgas / (gp->mmw * gp->density)); // otherwise we use the Jeans length
	    }
	    if (chemistry->H2_self_shielding >= 3) {
	    /* H2 self-shielding following Wolcott-Green & Haiman (2019),
 	       range of validity: T=100-8000 K, n<=1e7 cm^-3 */
	        const double N_H2 = cunits.dom * gp->rhoH2 * l_H2shield;
                const double tgas_touse = min(max(gp->tgas,100.),8000);
                const double ngas_touse = min(gp->density * cunits.dom / gp->mmw, 1.e7);
                const double aWG2019 = (0.8711 * log10(tgas_touse) - 1.928) *
                    exp(-0.2856 * log10(ngas_touse)) +
                    (-0.9639 * log10(tgas_touse) + 3.892);
                const double x = 2.e-15 * N_H2;
                const double b_doppler = 1.e-5 * sqrt(2. * kboltz * gp->tgas / mh);
                const double f_shield = 0.965 / pow(1. + x/b_doppler, aWG2019) +
                    0.035 * exp(-8.5e-4 * sqrt(1. + x)) / sqrt(1. + x);

	        my_rates->k31shield *= min(f_shield, 1.);
	    }
	    else if (chemistry->H2_self_shielding == 4) {  
	        // H+H2 self-shielding from Schauer+15 eq 8,9
	        const double NH_cgs = cunits.dom * gp->rhoH * l_H2shield;
	        const double xH = NH_cgs / 2.85e23;
	        const double fH_shield = pow(1.+xH,-1.62) * exp(-0.149*xH);
	        const double NH2_cgs = cunits.dom * gp->rhoH2 * l_H2shield;
    	        const double DH2_cgs = 1.e-5 * sqrt(2.*1.38e-16*gp->tgas / 3.346e-24);
    	        const double xH2 = NH2_cgs / 8.465e13;
    	        const double fH2_shield = 0.9379/pow(1.f+xH2/DH2_cgs,1.879) + 0.03465/pow(1.f+xH2,0.473) * exp(-2.293e-4*sqrt(1+xH2));

	        my_rates->k31shield *= min(fH_shield, 1.) * min(fH2_shield, 1.);
		//if (gp->verbose) printf("k31: %g %g %g %g %g %g %g\n",my_rates->k31shield, NH_cgs, NH2_cgs, gp->nH, gp->rhoH2/mh, fH_shield, fH2_shield);
	    }
	}
	if (chemistry->H2_self_shielding > 0 && chemistry->H2_custom_shielding == 1) {
	    // user specifies the H2 shielding factor directly
	    my_rates->k31shield *= gp->H2_custom_shielding_factor;
	}
}

static inline void init_temperature_interpolation(grackle_part_data *gp, chemistry_data *chemistry, interp_struct *interpolation, crackle_units cunits, chemistry_data_storage gr)
{
	/* mean mol weight and specific thermal energy */
	double mmw_neutral = gp->density / (gp->HI_density + gp->HII_density + 0.25 * (gp->HeI_density + gp->HeII_density + gp->HeIII_density) + 0.5 * (gp->H2I_density + gp->H2II_density) + gp->HM_density);
        gp->u_cmb = TCMB0 * (1.+cunits.redshift) / (chemistry->Gamma - 1.) / cunits.temp_to_K / mmw_neutral;

	/* interpolation table setup */
        interpolation->logT_lo = log(chemistry->TemperatureStart);
        interpolation->logT_hi = log(chemistry->TemperatureEnd);
        interpolation->dlogT = (interpolation->logT_hi - interpolation->logT_lo)/(chemistry->NumberOfTemperatureBins-1);
        interpolation->dlogT_inv = 1./interpolation->dlogT;

	/* cloudy table setup for this particle */
	/* Set up the density lookup; this doesn't change during the cooling step */
	interpolation->cloudy_par[0] = log10(gp->nH);
	interpolation->cloudy_dxinv[0] = 1.f * (gr.cloudy_metal.grid_dimension[0]-1) / (gr.cloudy_metal.grid_parameters[0][gr.cloudy_metal.grid_dimension[0]-1] - gr.cloudy_metal.grid_parameters[0][0]);  // Inverse grid spacing
	interpolation->cloudy_index[0] = (int) ((interpolation->cloudy_par[0] - gr.cloudy_metal.grid_parameters[0][0]) * interpolation->cloudy_dxinv[0]);
	//interpolation->cloudy_delta[0] = ((interpolation->cloudy_par[0] - gr.cloudy_metal.grid_parameters[0][interpolation->cloudy_index[0]]) * interpolation->cloudy_dxinv[0]);
	if (interpolation->cloudy_index[0] < 0) interpolation->cloudy_index[0] = 0;
	if (interpolation->cloudy_index[0] > gr.cloudy_metal.grid_dimension[0]-2) interpolation->cloudy_index[0] = gr.cloudy_metal.grid_dimension[0]-2;

	/* Set up the redshift lookup; the redshift table is not uniformly spaced so we search for the right index */
	interpolation->cloudy_par[1] = cunits.redshift;
	for (int i=1; i<gr.cloudy_metal.grid_dimension[1]-1; i++) {
	    interpolation->cloudy_index[1] = i-1;
	    if (interpolation->cloudy_par[1] < gr.cloudy_metal.grid_parameters[1][i]) break;
	}
	interpolation->cloudy_dxinv[1] = 1. / (gr.cloudy_metal.grid_parameters[1][interpolation->cloudy_index[1]+1] - gr.cloudy_metal.grid_parameters[1][interpolation->cloudy_index[1]]);
	//interpolation->cloudy_delta[1] = ((interpolation->cloudy_par[1] - gr.cloudy_metal.grid_parameters[1][interpolation->cloudy_index[1]]) * interpolation->cloudy_dxinv[1]);

	/* Accuracy value must be non-zero otherwise code will not do cooling */
	assert(chemistry->accuracy > 0.f && "CRACKLE ERROR: chemistry->accuracy input value must be greater than zero! (typically 0.1)");

	return;
}

static inline void setup_temperature_interpolation(double tgas, chemistry_data *chemistry, interp_struct *interpolation)
{
        /* Compute indexing for lookup tables */
        interpolation->logT = log(tgas);
        if (interpolation->logT > interpolation->logT_hi) interpolation->logT = interpolation->logT_hi;
        if (interpolation->logT < interpolation->logT_lo) interpolation->logT = interpolation->logT_lo;
        interpolation->index = (int)((interpolation->logT-interpolation->logT_lo) * interpolation->dlogT_inv);
        if (interpolation->index < 0) interpolation->index = 0;  
	if (interpolation->index > chemistry->NumberOfTemperatureBins-2) interpolation->index = chemistry->NumberOfTemperatureBins-2;

	/* Compute fraction within T bin */
        const double t1 = interpolation->logT_lo + interpolation->index * interpolation->dlogT;
        interpolation->binfrac = (interpolation->logT - t1) * interpolation->dlogT_inv;

	return;
}

static inline double interpolate_rates(double *r, double b, int i)
{
	return r[i] + b * (r[i+1]-r[i]);
}

static inline void primordial_cooling_rates(interp_struct *interpolation, chemistry_data_storage grackle_rates, chemistry_rate_storage *my_rates)
{
	my_rates->ceHI = interpolate_rates(grackle_rates.ceHI, interpolation->binfrac, interpolation->index);
	my_rates->ceHeI = interpolate_rates(grackle_rates.ceHeI, interpolation->binfrac, interpolation->index);
	//printf("ceHeI: %g %g %g %d %g\n",grackle_rates.ceHeI[193],grackle_rates.ceHeI[194],my_rates->ceHeI,interpolation->index,interpolation->binfrac);
	my_rates->ceHeII = interpolate_rates(grackle_rates.ceHeII, interpolation->binfrac, interpolation->index);
	my_rates->ciHI = interpolate_rates(grackle_rates.ciHI, interpolation->binfrac, interpolation->index);
	my_rates->ciHeI = interpolate_rates(grackle_rates.ciHeI, interpolation->binfrac, interpolation->index);
	my_rates->ciHeIS = interpolate_rates(grackle_rates.ciHeIS, interpolation->binfrac, interpolation->index);
	my_rates->ciHeII = interpolate_rates(grackle_rates.ciHeII, interpolation->binfrac, interpolation->index);
	my_rates->reHII = interpolate_rates(grackle_rates.reHII, interpolation->binfrac, interpolation->index);
	my_rates->reHeII1 = interpolate_rates(grackle_rates.reHeII1, interpolation->binfrac, interpolation->index);
	my_rates->reHeII2 = interpolate_rates(grackle_rates.reHeII2, interpolation->binfrac, interpolation->index);
	my_rates->reHeIII = interpolate_rates(grackle_rates.reHeIII, interpolation->binfrac, interpolation->index);
	my_rates->brem = interpolate_rates(grackle_rates.brem, interpolation->binfrac, interpolation->index);
}


static inline void H2_cooling_rates(interp_struct *interpolation, chemistry_data_storage grackle_rates, chemistry_rate_storage *my_rates)
{
	my_rates->GAHI = interpolate_rates(grackle_rates.GAHI, interpolation->binfrac, interpolation->index);
	my_rates->GAH2 = interpolate_rates(grackle_rates.GAH2, interpolation->binfrac, interpolation->index);
	my_rates->GAHe = interpolate_rates(grackle_rates.GAHe, interpolation->binfrac, interpolation->index);
	my_rates->GAHp = interpolate_rates(grackle_rates.GAHp, interpolation->binfrac, interpolation->index);
	my_rates->GAel = interpolate_rates(grackle_rates.GAel, interpolation->binfrac, interpolation->index);
	my_rates->H2LTE = interpolate_rates(grackle_rates.H2LTE, interpolation->binfrac, interpolation->index);
	my_rates->cieco = interpolate_rates(grackle_rates.cieco, interpolation->binfrac, interpolation->index);
}


static inline void primordial_species_rates(interp_struct *interpolation, chemistry_data_storage grackle_rates, chemistry_rate_storage *my_rates)
{
	my_rates->k1 = interpolate_rates(grackle_rates.k1, interpolation->binfrac, interpolation->index);
	my_rates->k2 = interpolate_rates(grackle_rates.k2, interpolation->binfrac, interpolation->index);
	my_rates->k3 = interpolate_rates(grackle_rates.k3, interpolation->binfrac, interpolation->index);
	my_rates->k4 = interpolate_rates(grackle_rates.k4, interpolation->binfrac, interpolation->index);
	my_rates->k5 = interpolate_rates(grackle_rates.k5, interpolation->binfrac, interpolation->index);
	my_rates->k6 = interpolate_rates(grackle_rates.k6, interpolation->binfrac, interpolation->index);
	my_rates->k57 = interpolate_rates(grackle_rates.k57, interpolation->binfrac, interpolation->index);
	my_rates->k58 = interpolate_rates(grackle_rates.k58, interpolation->binfrac, interpolation->index);
}

static inline void molecular_species_rates(interp_struct *interpolation, chemistry_data_storage grackle_rates, chemistry_rate_storage *my_rates)
{
	/* 9-species model for mode > 1 */
	my_rates->k7 = interpolate_rates(grackle_rates.k7, interpolation->binfrac, interpolation->index);
	my_rates->k8 = interpolate_rates(grackle_rates.k8, interpolation->binfrac, interpolation->index);
	my_rates->k9 = interpolate_rates(grackle_rates.k9, interpolation->binfrac, interpolation->index);
	my_rates->k10 = interpolate_rates(grackle_rates.k10, interpolation->binfrac, interpolation->index);
	my_rates->k11 = interpolate_rates(grackle_rates.k11, interpolation->binfrac, interpolation->index);
	my_rates->k12 = interpolate_rates(grackle_rates.k12, interpolation->binfrac, interpolation->index);
	my_rates->k13 = interpolate_rates(grackle_rates.k13, interpolation->binfrac, interpolation->index);
	my_rates->k14 = interpolate_rates(grackle_rates.k14, interpolation->binfrac, interpolation->index);
	my_rates->k15 = interpolate_rates(grackle_rates.k15, interpolation->binfrac, interpolation->index);
	my_rates->k16 = interpolate_rates(grackle_rates.k16, interpolation->binfrac, interpolation->index);
	my_rates->k17 = interpolate_rates(grackle_rates.k17, interpolation->binfrac, interpolation->index);
	my_rates->k18 = interpolate_rates(grackle_rates.k18, interpolation->binfrac, interpolation->index);
	my_rates->k19 = interpolate_rates(grackle_rates.k19, interpolation->binfrac, interpolation->index);
	my_rates->k22 = interpolate_rates(grackle_rates.k22, interpolation->binfrac, interpolation->index);
	/* H2 formation heating terms */
	my_rates->n_cr_n = interpolate_rates(grackle_rates.n_cr_n, interpolation->binfrac, interpolation->index);
	my_rates->n_cr_d1 = interpolate_rates(grackle_rates.n_cr_d1, interpolation->binfrac, interpolation->index);
	my_rates->n_cr_d2 = interpolate_rates(grackle_rates.n_cr_d2, interpolation->binfrac, interpolation->index);
}

static inline void deuterium_species_rates(int chemistry_flag, interp_struct *interpolation, chemistry_data_storage grackle_rates, chemistry_rate_storage *my_rates)
{
        if (chemistry_flag > 2) {
	    /* 12-species model for mode > 2 */
	    my_rates->k50 = interpolate_rates(grackle_rates.k50, interpolation->binfrac, interpolation->index);
	    my_rates->k51 = interpolate_rates(grackle_rates.k51, interpolation->binfrac, interpolation->index);
	    my_rates->k52 = interpolate_rates(grackle_rates.k52, interpolation->binfrac, interpolation->index);
	    my_rates->k53 = interpolate_rates(grackle_rates.k53, interpolation->binfrac, interpolation->index);
	    my_rates->k54 = interpolate_rates(grackle_rates.k54, interpolation->binfrac, interpolation->index);
	    my_rates->k55 = interpolate_rates(grackle_rates.k55, interpolation->binfrac, interpolation->index);
	    my_rates->k56 = interpolate_rates(grackle_rates.k56, interpolation->binfrac, interpolation->index);
	}
	else {
	    my_rates->k50 = 0.;
	    my_rates->k51 = 0.;
	    my_rates->k52 = 0.;
	    my_rates->k53 = 0.;
	    my_rates->k54 = 0.;
	    my_rates->k55 = 0.;
	    my_rates->k56 = 0.;
	}
}

static inline void dust_species_rates(double tdust, double dust2gas, chemistry_data *chemistry, chemistry_data_storage grackle_rates, chemistry_rate_storage *my_rates, interp_struct *interpolation)
{
	interp_struct dust_interp[1];

	/* limit temperature table bounds */
        dust_interp->logT = log(tdust);
        dust_interp->logT_lo = log(chemistry->DustTemperatureStart);
        dust_interp->logT_hi = log(chemistry->DustTemperatureEnd);

	/* interpolation table setup */
        dust_interp->dlogT = (dust_interp->logT_hi - dust_interp->logT_lo)/(chemistry->NumberOfDustTemperatureBins-1);
        dust_interp->dlogT_inv = 1./dust_interp->dlogT;
	if (dust_interp->logT < dust_interp->logT_lo + dust_interp->dlogT) dust_interp->logT = dust_interp->logT_lo + dust_interp->dlogT;
	if (dust_interp->logT > dust_interp->logT_hi) dust_interp->logT = dust_interp->logT_hi;

        /* Compute indexing for lookup tables */
        dust_interp->index = (int)((dust_interp->logT-dust_interp->logT_lo) * dust_interp->dlogT_inv);
        if (dust_interp->index < 1) dust_interp->index = 1;  
	if (dust_interp->index > chemistry->NumberOfDustTemperatureBins-1) dust_interp->index = chemistry->NumberOfDustTemperatureBins-1;

	/* Compute fraction within T bin */
        const double t1 = dust_interp->index * dust_interp->dlogT + dust_interp->logT_lo;
        dust_interp->binfrac = (dust_interp->logT - t1) * dust_interp->dlogT_inv;
	assert(dust_interp->binfrac >= 0.f);
	assert(dust_interp->binfrac <= 1.f);

	/* 2-D interpolation of dust-H2 formation rate table */
	const int index1 = interpolation->index + chemistry->NumberOfTemperatureBins * (dust_interp->index-1);
	const double dust1 = interpolate_rates(grackle_rates.h2dust, interpolation->binfrac, index1);
	const int index2 = interpolation->index + chemistry->NumberOfTemperatureBins * (dust_interp->index);
	const double dust2 = interpolate_rates(grackle_rates.h2dust, interpolation->binfrac, index2);

	my_rates->h2dust = (dust1 + (dust2 - dust1) * dust_interp->binfrac) * dust2gas;
	my_rates->gammah = grackle_rates.gammah;
	//printf("h2dust: %g %g %g\n",my_rates->h2dust, dust2gas, dust1);
}

static inline void lookup_chemistry_coeffs(int chemistry_flag, chemistry_data_storage grackle_rates, chemistry_rate_storage *my_rates, interp_struct *interpolation)
{
        /* Interpolate primordial species rate coefficients for given T */
        primordial_species_rates(interpolation, grackle_rates, my_rates);

        /* Interpolate molecular species rate coefficients for given T */
        if (chemistry_flag > 1) {
            molecular_species_rates(interpolation, grackle_rates, my_rates);
        }

        /* Interpolate deuterium species rate coefficients for given T */
        deuterium_species_rates(chemistry_flag, interpolation, grackle_rates, my_rates);
}

static inline double interpolate_3d(double *par, int *ind, double *cloudy_data, long long int *griddim, double *gridpar[3])
{
	int q,w,int_index;
	double value,value2[3], value3[3];
	double slope;
	for (q=1; q<3;q++) {
	    for (w=1; w<3;w++) {
		int_index = ((q+ind[0]-1)*griddim[1] + (w+ind[1]-1))*griddim[2] + (ind[2]);
		slope = (cloudy_data[int_index+1] - cloudy_data[int_index]) / (gridpar[2][ind[2]+1] - gridpar[2][ind[2]]);
		value3[w] = (par[2] - gridpar[2][ind[2]]) * slope + cloudy_data[int_index];
	    }
	    slope = (value3[2]-value3[1]) / log((1.+gridpar[1][ind[1]+1]) / (1.+gridpar[1][ind[1]]));
	    value2[q] = log((1.+par[1]) / (1.+gridpar[1][ind[1]])) * slope + value3[1];
	}
	slope = (value2[2]-value2[1]) / (gridpar[0][ind[0]+1] - gridpar[0][ind[0]]);
	value = (par[0] - gridpar[0][ind[0]]) * slope + value2[1];
	//assert( (par[2]>=0.9999*gridpar[2][ind[2]]) && (par[2]<=1.0001*gridpar[2][ind[2]+1]) ); // some slop to account for the cloudy table in floats vs crackle using doubles
	if (value != value) {
            printf("Something wrong in interpolate_3d! value=nan\n");
            pause();
        }
	return value;
}

static inline double cloudy_metal_cooling(double logtem, chemistry_data *chemistry, interp_struct *interpolation, chemistry_data_storage gr, crackle_units cunits)
{ 
	double log_edot, edot_metal = 0.f;

	/* variable translation to fortran code in cool_cloudy_g.F:
	 * clGridDim = gr.cloudy_metal.grid_dimension
	 * gridPar = gr.cloudy_metal.grid_parameters
	 * dgridPar = gr.cloudy_metal.grid_parameters
	 * input1,2,3 = par[3]
	 * dgridPar1,2,3 = dclPar(3) = 1/dgridinv[3]
	 * clDataSize = gr.cloudy_metal.grid_dimension[0]*gr.cloudy_metal.grid_dimension[1]*gr.cloudy_metal.grid_dimension{2]
	 * clCooling = gr.cloudy_metal.cooling_data
	 */

	/* The temperature changes each iternation, so we need to find the current index etc */
	interpolation->cloudy_par[2] = logtem * 0.4342944819; // constant to convert from log to log10
	if (interpolation->cloudy_par[2] < gr.cloudy_metal.grid_parameters[2][0]) interpolation->cloudy_par[2] = gr.cloudy_metal.grid_parameters[2][0];
	if (interpolation->cloudy_par[2] > gr.cloudy_metal.grid_parameters[2][gr.cloudy_metal.grid_dimension[2]-1]) interpolation->cloudy_par[2] = gr.cloudy_metal.grid_parameters[2][gr.cloudy_metal.grid_dimension[2]-1];
	interpolation->cloudy_dxinv[2] = 1.f * (gr.cloudy_metal.grid_dimension[2]-1) / (gr.cloudy_metal.grid_parameters[2][gr.cloudy_metal.grid_dimension[2]-1] - gr.cloudy_metal.grid_parameters[2][0]);  
	interpolation->cloudy_index[2] = (int) ((interpolation->cloudy_par[2] - gr.cloudy_metal.grid_parameters[2][0]) * interpolation->cloudy_dxinv[2]); /* Find grid index */
	if (interpolation->cloudy_index[2] < 0) interpolation->cloudy_index[2] = 0;
	if (interpolation->cloudy_index[2] > gr.cloudy_metal.grid_dimension[2]-2) interpolation->cloudy_index[2] = gr.cloudy_metal.grid_dimension[2]-2;
	//interpolation->cloudy_delta[2] = ((interpolation->cloudy_par[2] - gr.cloudy_metal.grid_parameters[2][interpolation->cloudy_index[2]]) * interpolation->cloudy_dxinv[2]); 

	/* Metal cooling */
	if (chemistry->metal_cooling) {
	    log_edot = interpolate_3d(interpolation->cloudy_par, interpolation->cloudy_index, gr.cloudy_metal.cooling_data, gr.cloudy_metal.grid_dimension, gr.cloudy_metal.grid_parameters);
	    edot_metal -= pow(10.f, log_edot);
	}

	/* CMB heating, if near CMB temperature */
	if (interpolation->cloudy_par[2] - cunits.log10_TCMB < 2.f) {
	    interpolation->cloudy_par[2] = cunits.log10_TCMB;  
	    interpolation->cloudy_dxinv[2] = 1.f * (gr.cloudy_metal.grid_dimension[2]-1) / (gr.cloudy_metal.grid_parameters[2][gr.cloudy_metal.grid_dimension[2]-1] - gr.cloudy_metal.grid_parameters[2][0]);  
	    interpolation->cloudy_index[2] = (int) ((interpolation->cloudy_par[2] - gr.cloudy_metal.grid_parameters[2][0]) * interpolation->cloudy_dxinv[2]); 
	    //interpolation->cloudy_delta[2] = ((interpolation->cloudy_par[2] - gr.cloudy_metal.grid_parameters[2][interpolation->cloudy_index[2]]) * interpolation->cloudy_dxinv[2]); 
	    if (interpolation->cloudy_index[2] < 0) interpolation->cloudy_index[2] = 0;
	    if (interpolation->cloudy_index[2] > gr.cloudy_metal.grid_dimension[2]-2) interpolation->cloudy_index[2] = gr.cloudy_metal.grid_dimension[2]-2;
	    log_edot = interpolate_3d(interpolation->cloudy_par, interpolation->cloudy_index, gr.cloudy_metal.cooling_data, gr.cloudy_metal.grid_dimension, gr.cloudy_metal.grid_parameters);
	    edot_metal += pow(10.f, log_edot);
	}

	/* UVB heating */
	if (chemistry->UVbackground) {
	    log_edot = interpolate_3d(interpolation->cloudy_par, interpolation->cloudy_index, gr.cloudy_metal.heating_data, gr.cloudy_metal.grid_dimension, gr.cloudy_metal.grid_parameters);
	    edot_metal += pow(10.f, log_edot);
	}

	if (edot_metal != edot_metal) {
                printf("Something wrong in cloudy_metal_cooling! edot_metal=nan\n");
                pause();
        }

	return edot_metal;
}

