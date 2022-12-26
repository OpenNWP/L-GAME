#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

lgame_home_dir=~/code/L-GAME
run_id=grid_generation
export OMP_NUM_THREADS=4

cat > namelist.nml << EOF

&run
run_id="$run_id"
lat_center=0.8929595951304794
lon_center=0.1199133716060684
ny=35
nx=35
dy=25000
dx=25000
run_span_min=0
/

&io
lread_oro=.false.
lwrite_grid=.true.
lread_land_sea=.false.
lset_oro=.true.
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
/

&bc
/

EOF

# That's it, here we go. Do not change anything below this line.
if [ ! -f $lgame_home_dir/grids/phys_sfc_properties/etopo.nc ]
then
  cd $lgame_home_dir/grids/phys_sfc_properties
  ./download_etopo.sh
  cd $lgame_home_dir
fi

source $lgame_home_dir/run_scripts/.sh/root_script.sh

cd $lgame_home_dir/grids/phys_sfc_properties
echo "Creating land fraction ..."
python3 land_fraction.py
echo "Land fraction mask created."

cd $lgame_home_dir













