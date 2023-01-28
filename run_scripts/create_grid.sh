#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

lgame_home_dir=~/code/L-GAME
run_id=grid_generation
export OMP_NUM_THREADS=4
orography_id=1

cat > namelist.nml << EOF

&run
run_id="$run_id"
lat_center=0.8929595951304794
lon_center=0.1199133716060684
ny=35
nx=37
dy=25000
dx=25000
run_span_min=0
/

&io
lread_geo=.false.
lwrite_grid=.true.
/

&constituents
lmoist=.false.
/

&diff
/

&rad
lrad=.false.
/

&surface
orography_id=$orography_id
lsleve=.true.
/

&bc
/

EOF

# That's it, now the basic run script will be sourced. Do not change anything below this line.

# creating needed directories
if [ ! -d $lgame_home_dir/grids ]
then
  mkdir $lgame_home_dir/grids
fi

if [ ! -d $lgame_home_dir/grids/phys_sfc_quantities ]
then
  mkdir $lgame_home_dir/grids/phys_sfc_quantities
fi

# downloading land use data if necessary
if [ $orography_id -eq 1 ]
then
  source $lgame_home_dir/run_scripts/.sh/download_phys_sfc_quantities.sh
fi

source $lgame_home_dir/run_scripts/.sh/root_script.sh











