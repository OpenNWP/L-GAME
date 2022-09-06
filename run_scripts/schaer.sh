#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

lgame_home_dir=~/code/L-GAME
run_id=schaer
export OMP_NUM_THREADS=4

cat > namelist.nml << EOF

&run
run_id="$run_id"
ny=3
nx=401
nlays=65
nlays_oro=65
toa=19500.0
dy=500
dx=500
run_span_min=360
scenario="schaer"
lcorio=.false.
sigma=1.0
lplane=.true.
lmoist=.false.
/

&io
lread_oro=.false.
lread_land_sea=.false.
dt_write_min=10
/

&constituents
lmoist=.false.
/

&diff
lmom_diff_h=.true.
diff_h_smag_rot=0.0
lmom_diff_v=.false.
ltemp_diff_h=.true.
ltemp_diff_v=.false.
lmass_diff_h=.false.
lmass_diff_v=.false.
/

&rad
lrad=.false.
/

&surface
orography_id=2
lprog_soil_temp=.false.
lsfc_sensible_heat_flux=.false.
lsfc_phase_trans=.false.
/

&bc
lperiodic=.true.
/

EOF

# That's it, here we go.
source $lgame_home_dir/run_scripts/.sh/root_script.sh












