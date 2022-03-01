import numpy as np
import netCDF4 as nc
import sys
from global_land_mask import globe

input_filename = "../grid.nc"
output_filename = "is_land.nc"

# end of the usual input section

# reading the model grid
ds = nc.Dataset(input_filename, "r", format = "NETCDF4")
lat_array = ds["lat_geo"][:, :]
lon_array = ds["lon_geo"][:, :]
ds.close()

is_land = np.zeros([len(lat_array[:, 0]), len(lat_array[0, :])], dtype=np.int8)

for i in range(len(is_land[:, 0])):
	for j in range(len(is_land[0, :])):
		if globe.is_land(np.rad2deg(lat_array[i, j]), np.rad2deg(lon_array[i, j])):
			is_land[i] = 1

output_filename = output_filename
ds = nc.Dataset(output_filename, "w", format = "NETCDF4")
ds.createDimension("lat", len(is_land[:, 0]))
ds.createDimension("lon", len(is_land[0, :]))
is_land_nc = ds.createVariable("is_land", int, ("lat", "lon"))
is_land_nc[:, :] = is_land
ds.close()

