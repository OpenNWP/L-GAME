# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

import sys
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc

run_id = sys.argv[1]
plot_time_since_init_min = sys.argv[2]
varname = sys.argv[3]

# reading the model output
input_filename = "../output/" + run_id + "/" + run_id + "+" + plot_time_since_init_min + "min.nc"
ds = nc.Dataset(input_filename, "r", format="NETCDF4")
lat_vector = ds["lat_model"][:]
lon_vector = ds["lon_model"][:]
z_array = ds["z"][:]
# reading the variable
plot_array = ds[varname][:]
unit = ds[varname].Unit
ds.close()

# preparations for the plot
semiminor = 6356752.314
semimajor = 6378137.0
re = (semimajor*semimajor*semiminor)**(1/3)
x_vector = lon_vector*re
x_array = np.zeros([len(z_array[:,0,1]), len(z_array[0,:,1])])
for i in range(len(x_array[:, 0])):
	x_array[i, :] = x_vector

# plotting
delta = np.max(plot_array[:,:,1]) - np.min(plot_array[:,:,1])
bounds=np.linspace(np.min(plot_array[:,:,1]) - 0.1*delta, np.max(plot_array[:,:,1]) + 0.1*delta, 1000)
fig = plt.figure()
cf = plt.contourf(1e-3*x_array[:,50:150], 1e-3*z_array[:,50:150,1], plot_array[:,50:150,1], cmap = "jet", levels = bounds)
plt.colorbar(cf, label = unit)
plt.title(run_id + " + " + plot_time_since_init_min + " min, var: " + varname)
plt.ylim([0, 10])
plt.xlabel("x / km")
plt.ylabel("z / km")
fig.savefig("../figs/" + run_id + "+" + plot_time_since_init_min + "min_" + varname + ".png")
