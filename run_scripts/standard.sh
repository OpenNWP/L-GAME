#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

lgame_home_dir=~/code/L-GAME
run_id=standard
export OMP_NUM_THREADS=2

cat > namelist.nml << EOF

&run
run_id="$run_id"
nlins=51
ncols=51
nlays=40
nlays_oro=30
dy=48000
dx=52000
run_span_hr=60
scenario="standard"
/

&diff
lklemp=.true.
lmom_diff_h=.false.
lmom_diff_v=.false.
/

&constituents
/

&surface
/

&bc
/

&io
dt_write_min=1
/

EOF

# That's it, here we go.
source $lgame_home_dir/run_scripts/.sh/root_script.sh
