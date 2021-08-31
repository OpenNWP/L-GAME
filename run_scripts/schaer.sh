#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

lgame_home_dir=~/code/L-GAME
run_id=schaer
export OMP_NUM_THREADS=2

cat > namelist.nml << EOF

&run
run_id="$run_id"
nlins=1
ncols=201
nlays=65
toa=19500.0
dy=1000
dx=1000
run_span_hr=10
scenario="schaer"
lcorio=.false.
l_z_equidist=.true.
/

&diff
lklemp=.true.
lmom_diff_h=.false.
lmom_diff_v=.false.
/

&io
dt_write_min=1
/

EOF

# That's it, here we go.
source .sh/root_script.sh
