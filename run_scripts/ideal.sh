#!/bin/bash

# This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/RFPET

rfpet_home_dir=~/code/rfpet
run_id=ideal
ncpus=1

cat > namelist.nml << EOF

&run
dy=500
dtime=100
dx=550
nlins=101
ncols=101
nlevs=100
run_span_hr=1
/

EOF

# That's it, here we go.
source root_script.sh
