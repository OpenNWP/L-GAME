#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

lgame_home_dir=~/code/L-GAME
run_id=advection
export OMP_NUM_THREADS=2

cat > namelist.nml << EOF

&run
run_id="$run_id"
nlins=3
ncols=301
nlays=50
nlays_oro=50
toa=25000.0
dy=1000
dx=1000
run_span_hr=2
scenario="advection"
lcorio=.false.
sigma=1.0
lplane=.true.
/

&diff
/

&constituents
no_of_gaseous_constituents=2
no_of_condensed_constituents=4
/

&surface
orography_id=3
lsoil_heat_conduction=.false.
lsfc_phase_trans=.false.
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
dt_write_min=60
/

EOF

# That's it, here we go.
source $lgame_home_dir/run_scripts/.sh/root_script.sh












