#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

lgame_home_dir=~/code/L-GAME
run_id=standard
export OMP_NUM_THREADS=4

cat > namelist.nml << EOF

&run
run_id="$run_id"
lat_center=0.8929595951304794
lon_center=0.1199133716060684
ny=35
nx=35
dy=25000
dx=25000
run_span_min=4320
/

&io
lread_geo=.true.
lwrite_grid=.false.
/

&constituents
lmoist=.true.
/

&diff
/

&rad
lrad=.true.
/

&surface
orography_id=1
lsleve=.true.
/

&bc
lperiodic=.true.
/

EOF

# That's it, here we go.
source $lgame_home_dir/run_scripts/.sh/root_script.sh














