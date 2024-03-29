# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

# This file is for plotting vertical slices.

import sys
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc

run_id = "schaer"
plot_time_since_init_min = 300
varname = "w"
unit_string = "m/s"
lgame_home_dir = "/home/max/code/L-GAME"

# end of usual input section

plot_time_since_init_min_str = str(plot_time_since_init_min)

# reading the model output
input_filename = lgame_home_dir + "/output/" + run_id + "/" + run_id + "+" + plot_time_since_init_min_str + "min.nc"
ds = nc.Dataset(input_filename, "r", format="NETCDF4")
lat_vector = ds["lat_model"][:]
lon_vector = ds["lon_model"][:]
z_array = ds["z"][:]
# reading the variable
plot_array = ds[varname][:]
unit = ds[varname].Unit
ds.close()

# preparations for the plot
x_vector = lon_vector*6371000.789927
x_array = np.zeros([len(z_array[:,0,1]), len(z_array[0,:,1])])
for i in range(len(x_array[:, 0])):
	x_array[i, :] = x_vector

# plotting
if run_id == "schaer":
	bounds = np.arange(-0.5, 0.5, 0.05)
	bounds = bounds[abs(bounds) > 0.01]
fig = plt.figure()
if run_id == "schaer":
	plt.contour(1e-3*x_array[:,150:250], 1e-3*z_array[:,150:250,1], plot_array[:,150:250,1], levels = bounds, colors = "black")
else:
	cf = plt.contourf(1e-3*x_array[:,150:250], 1e-3*z_array[:,150:250,1], plot_array[:,150:250,1])
	cbar = plt.colorbar(cf, fraction = 0.02, pad = 0.15, aspect = 80, orientation = "horizontal")
	cbar.ax.tick_params(labelsize = 12)
	cbar.set_label(unit_string, fontsize = 12)
plt.title(run_id + " + " + plot_time_since_init_min_str + " min, var: " + varname)
plt.ylim([0, 10])
plt.xlabel("x / km")
plt.ylabel("z / km")
fig.savefig(lgame_home_dir + "/output/" + run_id + "/" + run_id + "+" + plot_time_since_init_min_str + "min_" + varname + ".png")





