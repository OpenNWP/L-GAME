#!/bin/bash

# This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/RFPET

rfpet_home_dir=~/code/rfpet
run_id=ideal
ncpus=1

cat > namelist_$run_id << EOF

&run_ctl
run_id=$run_id
run_span_hr=1
/

EOF

# That's it, here we go.
source root_script.sh
