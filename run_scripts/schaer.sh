#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

lgame_home_dir=~/code/L-GAME
run_id=schaer
export OMP_NUM_THREADS=4

cat > namelist.nml << EOF

&run
run_id="$run_id"
ny=3
nx=401
n_layers=70
n_oro_layers=69
toa=21000.0
dy=500
dx=500
run_span_min=300
scenario="schaer"
lcorio=.false.
sigma=1.0
lplane=.true.
theta_adv_order=3
/

&io
lread_oro=.false.
lread_land_sea=.false.
lwrite_integrals=.true.
dt_write_min=10
/

&constituents
lmoist=.false.
/

&diff
lmom_diff_h=.false.
lmom_diff_v=.false.
ltemp_diff_h=.false.
ltemp_diff_v=.false.
lmass_diff_h=.false.
lmass_diff_v=.false.
klemp_begin_rel=0.5
/

&rad
lrad=.false.
/

&surface
orography_id=2
lsleve=.true.
lprog_soil_temp=.false.
lsfc_sensible_heat_flux=.false.
lsfc_phase_trans=.false.
lpbl=.false.
/

&bc
lperiodic=.true.
lfreeslip=.true.
/

EOF

# That's it, here we go.
source $lgame_home_dir/run_scripts/.sh/root_script.sh












