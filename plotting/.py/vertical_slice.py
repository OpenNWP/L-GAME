# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

import sys;
import matplotlib.pyplot as plt;
import numpy as np;
import netCDF4 as nc;

run_id = sys.argv[1];
plot_time_since_init_min = sys.argv[2];
varname = sys.argv[3];

input_filename = "../output/" + run_id + "/" + run_id + "+" + plot_time_since_init_min + "min.nc";
ds = nc.Dataset(input_filename, "r", format="NETCDF4");
x_vector = np.linspace(0,100*800,101);
z_vector = np.linspace(0,20,20);
# reading the variable
plot_array = ds[varname][:];
ds.close();

fig = plt.figure();
plt.contour(1e-3*x_vector, 1e-3*z_vector, plot_array[:,:,2]);
plt.title(run_id + " + " + plot_time_since_init_min + " min, var: " + varname);
plt.xlabel("x / km");
plt.ylabel("z / km");
fig.savefig("../figs/" + run_id + "+" + plot_time_since_init_min + "min_" + varname + ".png");
