**NOTE: crackle is retired.** This repository is frozen at its final version and will not be developed further. Crackle now lives inside the SWIFT simulation code as part of the Kiara model (https://github.com/kiara-project/SWIFT, under src/cooling/KIARA/crackle), where it runs directly on the SWIFT particles and no longer needs Grackle. All further development (including multi-species dust) is happening there.

Crackle is a fully C-based version of the Grackle cooling and chemistry code.  Information on Grackle is available here:
https://grackle.readthedocs.io/en/latest/index.html

Currently crackle ONLY works with Ewan Jones's kiara branch of Grackle, which includes the subgrid ISM dust+H2 model described here: 
https://ui.adsabs.harvard.edu/abs/2024MNRAS.535.1293J/abstract
Ewan's Grackle branch is here:
https://github.com/EwanBJones98/grackle
Don't forget to checkout the kiara branch.

To use Crackle, replace/add the codes in this repo to those in the above Grackle repo under src/clib (see below).  Then when putting it in your C code, set use_grackle=2.  The API is otherwise the same.

INSTALLATION INSTRUCTIONS

Clone crackle into your /src/clib/ directory of grackle.  It should then appear in the crackle subdirectory.

There are 6 files in grackle that crackle overwrites, in addition to adding many new ones.  From the /src/lib  directory, run the following commands to archive them in a separate directory, and then make symbolic links to (or copy over) all the crackle files:

mkdir grackle_orig

mv ./calculate_cooling_time.c grackle_orig

mv ./calculate_temperature.c grackle_orig

mv ./calculate_pressure.c grackle_orig

mv ./grackle_chemistry_data.h grackle_orig

mv ./grackle_types.h grackle_orig

mv ./solve_chemistry.c grackle_orig

ln -s crackle/* .

Then proceed with making grackle as usual.  
To update crackle, just git pull the changes into the crackle subdirectory and all should proceed as normal (if you've made symlinks).

TESTING

After building grackle (make in src/clib, no need to install), run the regression tests with:

crackle/tests/run_tests.sh

This evolves a gas parcel in several regimes (hot CGM, UVB-heated, self-shielded, warm/cold/hot ISM gas) with one crackle call per timestep, and checks the temperature against a converged reference using 100 substeps at higher accuracy.  It also checks specific fixes (H2 formation on dust, RT H2 dissociation, calculate_temperature before the first solve).  Add -g to also print the original Fortran grackle results for comparison.  The Cloudy data files are taken from grackle's input directory; set CRACKLE_DATA_DIR to use another location.
