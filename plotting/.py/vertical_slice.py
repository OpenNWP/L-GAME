# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

import sys;
import matplotlib.pyplot as plt;
import numpy as np;

run_id = sys.argv[1];
plot_time_since_init_min = sys.argv[2];
varname = sys.argv[3];

fig = plt.figure();
plt.title("varname");
fig.savefig("../figs/" + run_id + "+" + plot_time_since_init_min + "min_" + varname + ".png");
