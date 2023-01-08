#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

lgame_home_dir=~/code/L-GAME
run_id=${BASH_ARGV[9]}
export OMP_NUM_THREADS=${BASH_ARGV[20]}

cat > namelist.nml << EOF

&run
run_id="$run_id"
ny=${BASH_ARGV[8]}
nx=${BASH_ARGV[7]}
n_layers=${BASH_ARGV[6]}
n_oro_layers=${BASH_ARGV[5]}
dy=${BASH_ARGV[4]}
dx=${BASH_ARGV[3]}
run_span_min=${BASH_ARGV[2]}
lat_center=${BASH_ARGV[1]}
lon_center=${BASH_ARGV[0]}
start_year=${BASH_ARGV[10]}
start_month=${BASH_ARGV[11]}
start_day=${BASH_ARGV[12]}
start_hour=${BASH_ARGV[13]}
start_minute=${BASH_ARGV[14]}
/

&io
oro_raw_filename=${BASH_ARGV[15]}
restart_filename=${BASH_ARGV[17]}
/

&constituents
/

&diff
/

&rad
/

&surface
/

&bc
dtime_bc=${BASH_ARGV[18]}
bc_root_filename=${BASH_ARGV[19]}
/

EOF

# That's it, here we go.
source $lgame_home_dir/run_scripts/.sh/root_script.sh














