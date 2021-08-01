#!/bin/bash

# This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/RFPET

rfpet_home_dir=~/code/rfpet
run_id=ideal
ncpus=1

cat > namelist.nml << EOF

&run
run_span_hr=1
dtime=10
nlins=10
ncols=10
nlev=10
/

EOF

# That's it, here we go.
source root_script.sh
