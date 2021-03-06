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
scenario="standard"
/

&diff
/

&surface
/

&bc
lperiodic=.true.
/

&rad
/

&io
/

EOF

# That's it, here we go.
source $lgame_home_dir/run_scripts/.sh/root_script.sh














