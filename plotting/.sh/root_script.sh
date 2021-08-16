# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

if [ ! -d ../figs ]
then
mkdir ../figs
fi

python3 .py/vertical_slice.py $run_id $plot_time_since_init_min $varname
