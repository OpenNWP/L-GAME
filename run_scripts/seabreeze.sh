#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

lgame_home_dir=~/code/L-GAME
run_id=seabreeze
export OMP_NUM_THREADS=2

cat > namelist.nml << EOF

&run
run_id="$run_id"
start_year=2022
start_month=6
start_day=30
start_hour=0
start_minute=0
lat_center=0.7853981633974483
lon_center=0.0
nlins=3
ncols=101
nlays=50
nlays_oro=40
toa=20000.0
dy=2000
dx=2000
run_span_hr=24
scenario="seabreeze"
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
orography_id=0
/

&bc
lperiodic=.true.
/

&rad
dtime_rad=1800.0
/

&io
lread_oro=.false.
lread_land_sea=.false.
dt_write_min=60
/

EOF

# That's it, here we go.
source $lgame_home_dir/run_scripts/.sh/root_script.sh












