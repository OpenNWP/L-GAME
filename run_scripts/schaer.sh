#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

lgame_home_dir=~/code/L-GAME
run_id=schaer
export OMP_NUM_THREADS=2

cat > namelist.nml << EOF

&run
run_id="$run_id"
nlins=1
ncols=401
nlays=65
nlays_oro=65
toa=19500.0
dy=500
dx=500
run_span_hr=6
scenario="schaer"
lcorio=.false.
sigma=1.0
lplane=.true.
/

&diff
lmom_diff_h=.false.
lmom_diff_v=.false.
ltemp_diff_h=.false.
ltemp_diff_v=.false.
ltracer_diff_h=.false.
ltracer_diff_v=.false.
/

&constituents
no_of_gaseous_constituents=1
no_of_condensed_constituents=0
/

&surface
orography_id=2
lsoil=.false.
/

&bc
lperiodic=.true.
/

&rad
lrad=.false.
/

&io
lread_oro=.false.
lread_land_sea=.false.
dt_write_min=10
/

EOF

# That's it, here we go.
source $lgame_home_dir/run_scripts/.sh/root_script.sh












