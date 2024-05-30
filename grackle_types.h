/***********************************************************************
/
/ Grackle variable types
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef __GRACKLE_TYPES_H__
#define __GRACKLE_TYPES_H__
/***********************************************************************
/  
/ VARIABLE TYPES
/
************************************************************************/

#include "grackle_float.h"

#ifdef GRACKLE_FLOAT_4
#define gr_float float
#endif

#ifdef GRACKLE_FLOAT_8
#define gr_float double
#endif

/* Dust evolution definitions */
#define NUM_METAL_SPECIES_GRACKLE 10
#if defined(GRACKLE_FLOAT_4) == defined(GRACKLE_FLOAT_8)
#error "Both GRACKLE_FLOAT_4 and GRACKLE_FLOAT_8 are defined. Only one can be defined."
#endif

typedef struct
{

  int grid_rank;
  int *grid_dimension;
  int *grid_start;
  int *grid_end;

  gr_float grid_dx;

  gr_float *density;
  gr_float *HI_density;
  gr_float *HII_density;
  gr_float *HM_density;
  gr_float *HeI_density;
  gr_float *HeII_density;
  gr_float *HeIII_density;
  gr_float *H2I_density;
  gr_float *H2II_density;
  gr_float *DI_density;
  gr_float *DII_density;
  gr_float *HDI_density;
  gr_float *e_density;
  gr_float *metal_density;
  gr_float *dust_density;

  gr_float *internal_energy;
  gr_float *x_velocity;
  gr_float *y_velocity;
  gr_float *z_velocity;

  gr_float *volumetric_heating_rate;
  gr_float *specific_heating_rate;

  gr_float *temperature_floor;

  gr_float *RT_heating_rate;
  gr_float *RT_HI_ionization_rate;
  gr_float *RT_HeI_ionization_rate;
  gr_float *RT_HeII_ionization_rate;
  gr_float *RT_H2_dissociation_rate;

  gr_float *H2_self_shielding_length;
  gr_float *H2_custom_shielding_factor;

  gr_float *isrf_habing;
  
  // densities of individual metal species in gas phase
  gr_float *He_gas_metalDensity;
  gr_float *C_gas_metalDensity;
  gr_float *N_gas_metalDensity;
  gr_float *O_gas_metalDensity;
  gr_float *Ne_gas_metalDensity;
  gr_float *Mg_gas_metalDensity;
  gr_float *Si_gas_metalDensity;
  gr_float *S_gas_metalDensity;
  gr_float *Ca_gas_metalDensity;
  gr_float *Fe_gas_metalDensity;

  // densities of individual metal species in dust grains
  gr_float *He_dust_metalDensity;
  gr_float *C_dust_metalDensity;
  gr_float *N_dust_metalDensity;
  gr_float *O_dust_metalDensity;
  gr_float *Ne_dust_metalDensity;
  gr_float *Mg_dust_metalDensity;
  gr_float *Si_dust_metalDensity;
  gr_float *S_dust_metalDensity;
  gr_float *Ca_dust_metalDensity;
  gr_float *Fe_dust_metalDensity;
  
  gr_float *SNe_ThisTimeStep;

} grackle_field_data;

typedef struct
{

  int comoving_coordinates;
  double density_units;
  double length_units;
  double time_units;
  double velocity_units;
  double a_units;
  double a_value;
  double temperature_units;

} code_units;

typedef struct
{

  const char* version;
  const char* branch;
  const char* revision;

} grackle_version;

#endif