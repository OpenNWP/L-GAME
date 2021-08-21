#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

lgame_home_dir=~/code/L-GAME
run_id=schaer
ncpus=1

cat > namelist.nml << EOF

&run
run_id="$run_id"
nlins=3
ncols=201
nlays=65
toa=19500.0
dy=1050
dx=1000
run_span_hr=1
scenario="schaer"
llinear=.true.
lcorio=.false.
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

# That's it, here we go.
source root_script.sh