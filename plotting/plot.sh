#!/bin/bash

# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

lgame_home_dir=~/code/L-GAME # the home directory of L-GAME
run_id=schaer # the run ID of the run you want to plot
plot_time_since_init_min=300 # the time for which you want to have a plot
varname="w" # the variable you want to plot

source $lgame_home_dir/plotting/.sh/root_script.sh
