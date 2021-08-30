#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

lgame_home_dir=~/code/L-GAME
run_id=ideal
ncpus=1

cat > namelist.nml << EOF

&run
run_id="$run_id"
nlins=101
ncols=101
nlays=30
dy=850
dx=800
run_span_hr=1
scenario="standard"
/

&diff
lklemp=.false.
lmom_diff_h=.false.
lmom_diff_v=.false.
/

&io
dt_write_min=1
/

EOF

export OMP_NUM_THREADS=2

# That's it, here we go.
source root_script.sh
