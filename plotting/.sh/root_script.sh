# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

if [ ! -d $lgame_home_dir/figs ]
then
mkdir $lgame_home_dir/figs
fi

python3 $lgame_home_dir/plotting/.py/vertical_slice.py $run_id $plot_time_since_init_min $varname $lgame_home_dir
