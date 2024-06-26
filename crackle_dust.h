
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DUST_GROWTH_CAP 0.9 // Total dust growth is limited to this fraction of gas-phase metals
#define DUST_CHANGE_CAP 10.0  // factor by which dust density can go up or down in a single step
#define T_SPUTTERING 2.e4  // Temeprature below which sputtering is negligible

static inline void evolve_dust(grackle_part_data *gp, chemistry_data *chemistry, code_units *units, int ism_flag, double dtit)
{
        // dust growth via accretion and destruction via shock and thermal sputtering
        int k;
        double drho[NUM_METAL_SPECIES_GRACKLE];  // change in dust metal density
	double rhomet[NUM_METAL_SPECIES_GRACKLE];  // total metallicity including both gas and dust
	double drhos;  // dust mass loss via shock and sputtering
        double tau_ref=1.e20, tau_accr=1.e20, tau_accr0=1.e20, f_shocked=0., tau_sput=1.e20;
	double sec_per_year = 3.1536e7;

	if (gp->dust2gas <= tiny || gp->dust_density <= tiny) return;
        double dens_cgs = gp->density * units->density_units / (units->a_value*units->a_value*units->a_value);  // to physical cgs
        for (k=0; k<NUM_METAL_SPECIES_GRACKLE; k++) drho[k] = 0.;

	// Dust growth: metals accreted from gas to dust; only allowed when in ISM
	if (ism_flag) {
            tau_ref = chemistry->dust_growth_tauref * 1.e9 * sec_per_year / units->time_units;
            tau_accr0 = tau_ref * (chemistry->dust_growth_densref/dens_cgs) * sqrt(gp->tdust/gp->tgas); 
            for (k=0; k<NUM_METAL_SPECIES_GRACKLE; k++){
		gp->dust_metalDensity[k] = fmax(gp->dust_metalDensity[k], tiny);
                rhomet[k] = gp->gas_metalDensity[k] + gp->dust_metalDensity[k];
                tau_accr = tau_accr0 * gp->metallicity;
		if (rhomet[k] < tiny || tau_accr < tiny) continue;
                drho[k] = fmin((gp->gas_metalDensity[k] / rhomet[k]) * (gp->dust_metalDensity[k]/tau_accr) * dtit, DUST_GROWTH_CAP * gp->gas_metalDensity[k]); // growth in a single step capped at fraction of the available gas metals
            }
	}

	/* Dust destruction */
	if (gp->dust_density > tiny) {  // can't destroy if there is no dust
	    /* SNe shocks */
            if (gp->SNe_density <= tiny) drhos = 0.0;  // Total change in density from all destruction
            else{
	        const double mass_unit = units->density_units * units->length_units * units->length_units * units->length_units;
                const double rhogas = gp->density - gp->dust_density;
                const double Ms100 = 6800.0*chemistry->sne_coeff*(100.0/chemistry->sne_shockspeed)*(100.0/chemistry->sne_shockspeed) * SolarMass / mass_unit;  //code units, gas mass shocked per SNe (Sedov-Taylor phase)
                f_shocked = fmin( (Ms100 * gp->SNe_density) / (rhogas+tiny), 1.);  // fraction of mass shock-heated
                drhos = gp->dust_density * f_shocked * chemistry->dust_destruction_eff;  // some fraction of dust is destroyed in shocked gas
	if (gp->verbose) printf("dust: %g %g %g %g %g %g\n",gp->density, gp->nH, gp->SNe_density, rhogas, f_shocked, drhos);
	    }
	    /* sputtering */
	    if (gp->tgas > T_SPUTTERING) {  // negligibly small below this temperature
	        tau_sput = 1.7e8 * sec_per_year / units->time_units
                           * (chemistry->dust_grainsize / 0.1) * (1.e-27 / dens_cgs)
                           * (pow(2.e6 / gp->tgas, 2.5) + 1.0); // sputtering timescale, Tsai & Mathews (1995)
                //if (tau_sput > 0.) drhos += gp->dust_density / tau_sput * 3.0 * dtit;
                drhos += gp->dust_density * (1. - exp(-3.* fmin(dtit / tau_sput, 5.)));
	    }
	    /* Destruction is applied to all elements, in proportion to their metallicity */
	    for (k = 0; k < NUM_METAL_SPECIES_GRACKLE; k++) {
                drho[k] -= gp->dust_metalDensity[k] / gp->dust_density * drhos;
	    }
	}

	/* Now evolve the dust metal density */
	gp->dust_density = 0.;
	double old_dustDensity;
        for (k=0; k<NUM_METAL_SPECIES_GRACKLE; k++){
	    /* Limit change to not allow negative dust */
	    old_dustDensity = gp->dust_metalDensity[k];
	    gp->dust_metalDensity[k] += drho[k];
	    if (gp->dust_metalDensity[k] > DUST_CHANGE_CAP*old_dustDensity) {
	    	gp->dust_metalDensity[k] = DUST_CHANGE_CAP*old_dustDensity;
		drho[k] = gp->dust_metalDensity[k] - old_dustDensity;
	    }
	    if (gp->dust_metalDensity[k] < old_dustDensity/DUST_CHANGE_CAP) {
	    	gp->dust_metalDensity[k] = old_dustDensity/DUST_CHANGE_CAP;
		drho[k] = gp->dust_metalDensity[k] - old_dustDensity;
	    }
	    /* Move metals from gas to dust (or vice versa) */
	    gp->gas_metalDensity[k] -= drho[k];
	    gp->dust_density += gp->dust_metalDensity[k];
        }
	if (gp->verbose) printf("dust: tau_accr=%g tau_sput=%g snerho=%g drhos=%g rhog=%g td=%g tg=%g dust=%g\n",tau_accr0,tau_sput/3.,gp->SNe_density,drhos,gp->density,gp->tdust,gp->tgas,gp->dust_density);
	assert(gp->dust_density == gp->dust_density);

	return;
}

static inline double dust_thermal_balance(double tdust, double tgas, double trad, double trad4, double gamma_isrf, double gasgr, double nh)
{
        const double kgr1 = 4.e-4;
        const double kgr200 = 16.f; /* dust-rad coupling at sublimation temp */
        const double tsubl = 1500.f;  /* dust sublimation temperature */
        const double sigma_sb = 5.670373e-5;  /* Stefan-Boltzmann const, in erg/s/cm^2/K^4  */
        double kgr, sol, tdust4;

        /* Initialize: If input tdust<=0, return initial guess for tdust */
        if (tdust <= 0.f) {
            tdust = fmax(trad, pow(0.25 * gamma_isrf / sigma_sb / kgr1, 0.17f));
            return tdust;
        }

        /* Otherwise compute dust thermal balance, return heating - cooling rate */
        /* Compute dust-radiation coupling constant */
        if (tdust < 200.f) kgr = kgr1 * tdust * tdust;
        else if (tdust > tsubl) kgr = kgr200;
        else kgr = fmax(tiny, kgr200 * pow(tdust/tsubl, 12.f));
        /* Compute dust heating-cooling rate */
        tdust4 = tdust * tdust * tdust * tdust;
        sol = gamma_isrf + 4.f * sigma_sb * kgr * (trad4 - tdust4) + gasgr * nh * (tgas-tdust);
        //printf("sol: %g %g %g %g %g %g\n",sol, tdust, gamma_isrf, kgr, trad4-tdust4,tgas-tdust);
}

static inline double calculate_dust_temp(double tgas, double nh, double gasgr, double gamma_isrf, double trad, double td)
{
        /* Solve for Tdust from ISRF heating, CMB heating, and gas-grain cooling */
        const double tol = 1.e-3;
        double eps=1.e-2, tdold = -1.e20, dsol, slope; // for Newton-Raphson
        //double tdlo=trad, tdhi=tgas; // bisection limits
        const double trad4 = trad * trad * trad * trad;
        int iter = 0, maxiter = 100, sol;

        /* Solve for Tdust via Newton-Raphson */
        if (td<=0.f) td = dust_thermal_balance(td, tgas, trad, trad4, gamma_isrf, gasgr, nh); // initial guess
        while (fabs(td-tdold) > tol * td) {
            sol = dust_thermal_balance(td, tgas, trad, trad4, gamma_isrf, gasgr, nh);
            dsol = dust_thermal_balance((1.f+eps)*td, tgas, trad, trad4, gamma_isrf, gasgr, nh) - sol;
            tdold = td;
	    if (dsol > 0.) td = td - td * eps * sol / dsol;  // Newton-Raphson guess
	    else break;  // converged
            if (td < trad) td = trad;  // Limit tdust to [tCMB, tgas]
            if (td > tgas) td = tgas;
            if (sol * (sol+dsol) < 0.f) eps *= 0.5f;  // we have passed minimum; reduce eps
            if (++iter > maxiter) {
                printf("Crackle: Non-convergence in calculate_dust_temp(), returning tdust=%g (tdold=%g, dsol=%g eps=%g slope=%g)\n",td, tdold, dsol, eps, slope);
                break;
            }
        }

        return td;
}

