#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

lgame_home_dir=~/code/L-GAME
run_id=advection
export OMP_NUM_THREADS=4

cat > namelist.nml << EOF

&run
run_id="$run_id"
ny=3
nx=301
nlays=50
nlays_oro=50
toa=25000.0
dy=1000
dx=1000
run_span_min=120
scenario="advection"
lcorio=.false.
sigma=1.0
lplane=.true.
/

&io
lread_oro=.false.
lread_land_sea=.false.
dt_write_min=60
/

&constituents
/

&diff
/

&rad
lrad=.false.
/

&surface
orography_id=3
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












