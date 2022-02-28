#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

lgame_home_dir=~/code/L-GAME
run_id=standard
export OMP_NUM_THREADS=2

cat > namelist.nml << EOF

&run
run_id="$run_id"
lat_center_deg=51.16281607668721
lon_center_deg=6.870530100211603
x_dir_deg=90
nlins=25
ncols=25
nlays=50
nlays_oro=40
dy=500
dx=500
run_span_hr=0
scenario="create_grid"
/

&diff
/

&constituents
/

&surface
/

&bc
/

&rad
/

&io
lwrite_grid=.true.
/

EOF

# That's it, here we go.
source $lgame_home_dir/run_scripts/.sh/root_script.sh
