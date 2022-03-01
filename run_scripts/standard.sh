#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

lgame_home_dir=~/code/L-GAME
run_id=standard
export OMP_NUM_THREADS=2

cat > namelist.nml << EOF

&run
run_id="$run_id"
nlins=25
ncols=25
nlays=50
nlays_oro=40
dy=500
dx=500
run_span_hr=30
scenario="standard"
/

&diff
/

&constituents
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
