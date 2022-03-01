#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

lgame_home_dir=~/code/L-GAME
run_id=${BASH_ARGV[9]}
export OMP_NUM_THREADS=2

cat > namelist.nml << EOF

&run
run_id="$run_id"
nlins=${BASH_ARGV[8]}
ncols=${BASH_ARGV[7]}
nlays=${BASH_ARGV[6]}
nlays_oro=${BASH_ARGV[5]}
dy=${BASH_ARGV[4]}
dx=${BASH_ARGV[3]}
run_span_hr=${BASH_ARGV[2]}
lat_center=${BASH_ARGV[1]}
lon_center=${BASH_ARGV[0]}
start_year=${BASH_ARGV[10]}
start_month=${BASH_ARGV[11]}
start_day=${BASH_ARGV[12]}
start_hour=${BASH_ARGV[13]}
start_minute=${BASH_ARGV[14]}
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
dt_write_min=15
/

EOF

# That's it, here we go.
source $lgame_home_dir/run_scripts/.sh/root_script.sh
