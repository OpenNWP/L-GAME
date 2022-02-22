#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

lgame_home_dir=~/code/L-GAME
run_id=${BASH_ARGV[10]}
export OMP_NUM_THREADS=2

cat > namelist.nml << EOF

&run
run_id="$run_id"
nlins=${BASH_ARGV[9]}
ncols=${BASH_ARGV[8]}
nlays=${BASH_ARGV[7]}
nlays_oro=${BASH_ARGV[6]}
dy=${BASH_ARGV[5]}
dx=${BASH_ARGV[4]}
run_span_hr=${BASH_ARGV[3]}
scenario="standard"
lat_center_deg=${BASH_ARGV[2]}
lon_center_deg=${BASH_ARGV[1]}
x_dir_deg=${BASH_ARGV[0]}
/

&diff
lklemp=.true.
lmom_diff_h=.true.
lmom_diff_v=.true.
/

&constituents
/

&surface
/

&io
dt_write_min=60
/

EOF

# That's it, here we go.
source $lgame_home_dir/run_scripts/.sh/root_script.sh
