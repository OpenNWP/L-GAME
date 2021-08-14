#!/bin/bash

# This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/RFPET

rfpet_home_dir=~/code/RFPET
run_id=ideal
ncpus=1

cat > namelist.nml << EOF

&run
nlins=3
ncols=101
nlays=20
dy=800
dx=850
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

# That's it, here we go.
source root_script.sh
