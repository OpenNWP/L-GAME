#!/bin/bash

# This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/RFPET

rfpet_home_dir=~/code/rfpet
run_id=ideal
ncpus=1

cat > namelist.nml << EOF

&run
dy=800 ! at sea level
dx=850 ! at equator and sea level
nlins=101
ncols=101
nlays=80
run_span_hr=1
/

EOF

# That's it, here we go.
source root_script.sh
