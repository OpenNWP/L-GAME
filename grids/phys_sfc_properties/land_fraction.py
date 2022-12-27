import numpy as np
import netCDF4 as nc
import sys
from global_land_mask import globe

input_filename = "../grid.nc"
output_filename = "land_fraction.nc"

# end of the usual input section

# reading the model grid
ds = nc.Dataset(input_filename, "r", format = "NETCDF4")
lat_array = ds["lat_geo"][:, :]
lon_array = ds["lon_geo"][:, :]
lat_array_u = ds["lat_geo_u"][:, :]
lon_array_u = ds["lon_geo_u"][:, :]
lat_array_v = ds["lat_geo_v"][:, :]
lon_array_v = ds["lon_geo_v"][:, :]
ds.close()

land_fraction = np.zeros([len(lat_array[:, 0]), len(lat_array[0, :])])

for i in range(len(lat_array[:, 0])):
	for j in range(len(lat_array[0, :])):
		land_fraction_local = 0.0
		if globe.is_land(np.rad2deg(lat_array[j, i]), np.rad2deg(lon_array[j, i])):
			land_fraction_local = land_fraction_local + 1.0/5.0
		if globe.is_land(np.rad2deg(lat_array_u[j, i]), np.rad2deg(lon_array_u[j, i])):
			land_fraction_local = land_fraction_local + 1.0/5.0
		if globe.is_land(np.rad2deg(lat_array_u[j + 1, i]), np.rad2deg(lon_array_u[j + 1, i])):
			land_fraction_local = land_fraction_local + 1.0/5.0
		if globe.is_land(np.rad2deg(lat_array_v[j, i]), np.rad2deg(lon_array_v[j, i])):
			land_fraction_local = land_fraction_local + 1.0/5.0
		if globe.is_land(np.rad2deg(lat_array_v[j, i + 1]), np.rad2deg(lon_array_v[j, i + 1])):
			land_fraction_local = land_fraction_local + 1.0/5.0
		land_fraction[j, i] = land_fraction_local

output_filename = output_filename
ds = nc.Dataset(output_filename, "w", format = "NETCDF4")
ds.createDimension("lat", len(land_fraction[:, 0]))
ds.createDimension("lon", len(land_fraction[0, :]))
land_fraction_nc = ds.createVariable("land_fraction", float, ("lat", "lon"))
land_fraction_nc[:, :] = land_fraction
ds.close()

