""" functions to grid and re-grid netcdf files """ 

import xarray as xr
from scipy.interpolate import griddata
from netCDF4 import Dataset
import numpy as np


def interp_onto_regular_lat_lon_grid(da_in, grid_resolution_degrees, lat_var, lon_var, method='nearest', cheat=False):
    lons = np.array(da_in[lon_var].values)
    lats = np.array(da_in[lat_var].values)
    data = np.array(da_in.values)
    lons = lons.flatten()
    lats = lats.flatten()
    data = data.flatten()

    #Add wrap around at both longitude limits
    pts=np.squeeze(np.where(lons < -150))
    lons=np.append(lons, lons[pts]+360)
    lats=np.append(lats, lats[pts])
    data=np.append(data, data[pts])

    pts=np.squeeze(np.where(lons > 150))
    lons=np.append(lons, lons[pts]-360)
    lats=np.append(lats, lats[pts])
    data=np.append(data, data[pts])

    # Make the new grid
    lons_1d=np.arange(-180, 180, grid_resolution_degrees)
    lats_1d=np.arange(-90, 90, grid_resolution_degrees)
    lons_2d, lats_2d = np.meshgrid(lons_1d, lats_1d)
    
    # Interpolate Orca2 data onto the new regular grid
    data_new = griddata((lons, lats), data, (lons_2d, lats_2d), method=method)
    da_out = xr.DataArray(data=data_new, coords=[('lat', lats_1d),('lon', lons_1d)])
    ds_out = da_out.to_dataset(name=da_in.name)
    
    if cheat == False:
        ds_out[da_in.name].attrs = {'units':da_in.units, 'long_name':da_in.long_name}
    
    return ds_out

