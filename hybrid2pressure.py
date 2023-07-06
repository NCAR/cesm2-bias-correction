#This code mirrors the ncl script convert_cesm_hybric_nc_to_pressure_int.ncl
#The old code interpolates the grid, convert to pressure coordinates, and converts to WRF-INT format
#This code will only interpolate the grid and convert to pressure coordinates (no conversion to WRF-INT)

#INPUTS:    3D CESM 6hr data, in netCDF4 format
#OUTPUTS:   3D CESM 6hr data, interpolated onto lat/lon grid and in pressure coordinates, in netCDF4 format. 

#Author: Kenton Wu 
#Date: Jun 26, 2023 
#Email: wukenton@ucar.edu, wukenton@gmail.com, wukenton@utexas.edu 

import argparse 
import ast

import cf_xarray
import cftime
import geocat
import holoviews as hv
import hvplot
import hvplot.xarray
import intake
import numpy as np
import pop_tools
import xarray as xr
import xesmf as xe
from distributed import Client
from ncar_jobqueue import NCARCluster
from pop_tools.grid import _compute_corners


#Command line option handling ----------------------------------------------------------------------------------
parser = argparse.ArgumentParser()  

parser.add_argument('CASE',type=str, help='One of the following IPCC Climate Scenarios: 20THC/RCP85/RCP60/RCP45')
parser.add_argument('--o',type=str,help='Output directory path')

#File Handling ----------------------------------------------------------------------------------

print("Opening files...")

in_ta = xr.open_dataset("atmos_ta.nc")         # 6-hourly 3-d T
in_ua = xr.open_dataset("atmos_ua.nc")         # 6-hourly 3-d U
in_va = xr.open_dataset("atmos_va.nc")         # 6-hourly 3-d V
in_hus = xr.open_dataset("atmos_hus.nc")       # 6-hourly 3-d Q
in_ps = xr.open_dataset("atmos_ps.nc")         # 6-hourly surface pressure
in_zsfc = xr.open_dataset("atmos_zsfc.nc")     # static surface geopotential
in_lmask = xr.open_dataset("atmos_lmask.nc")   # static land mask
in_snw = xr.open_dataset("atmos_snw_1.nc")     # monthly SWE
in_mrlsl = xr.open_dataset("atmos_mrlsl_1.nc") # monthly soil moisture
in_ts = xr.open_dataset("atmos_ts_1.nc")       # monthly skin temp
in_tsl = xr.open_dataset("atmos_tsl_1.nc")     # monthly soil temp
in_tos = xr.open_dataset("atmos_tos_1.nc")     # daily SST on pop grid (gaussian)
in_sic = xr.open_dataset("atmos_sic_1.nc")     # daily SEAICE % on POP grid (gaussian)

#Interpolate SST and SEA ICE fields to CESM Atmospheric Domain ----------------------------------------------------------------------------------

print('Converting POP info into rectangular grid...')

#lat_corners, lon_corners = gen_corner_calc(era_data)

#Regrids SST grid to whatever the atmospheric grid is automatically
regrid = xe.Regridder(in_tos, in_ta, method = 'bilinear')
regrid.to_netcdf('gx1v6_to_era_latlon.nc')

rectangular_SST = regrid(in_tos)
print(rectangular_SST)
rectangular_SST.to_netcdf('python_regrid.nc')
#
#use some sort of broadcasting or view here to clone to a 6-hrly variable


#replace the default values here 
#def gen_corner_calc(ds, cell_corner_lat='ULAT', cell_corner_lon='ULONG'):
#    """
#    Generates corner information and creates single dataset with output
#    """
#
#    cell_corner_lat = ds[cell_corner_lat]
#    cell_corner_lon = ds[cell_corner_lon]
#    # Use the function in pop-tools to get the grid corner information
#    corn_lat, corn_lon = _compute_corners(cell_corner_lat, cell_corner_lon)
#
#    # Make sure this returns four corner points
#    assert corn_lon.shape[-1] == 4
#
#    lon_shape, lat_shape = corn_lon[:, :, 0].shape
#    out_shape = (lon_shape + 1, lat_shape + 1)
#
#    # Generate numpy arrays to store destination lats/lons
#    out_lons = np.zeros(out_shape)
#    out_lats = np.zeros(out_shape)
#
#    # Assign the northeast corner information
#    out_lons[1:, 1:] = corn_lon[:, :, 0]
#    out_lats[1:, 1:] = corn_lat[:, :, 0]
#
#    # Assign the northwest corner information
#    out_lons[1:, :-1] = corn_lon[:, :, 1]
#    out_lats[1:, :-1] = corn_lat[:, :, 1]
#
#    # Assign the southwest corner information
#    out_lons[:-1, :-1] = corn_lon[:, :, 2]
#    out_lats[:-1, :-1] = corn_lat[:, :, 2]
#
#    # Assign the southeast corner information
#    out_lons[:-1, 1:] = corn_lon[:, :, 3]
#    out_lats[:-1, 1:] = corn_lat[:, :, 3]
#
#    return out_lats, out_lons

#Turn monthly data into 6hr data ----------------------------------------------------------------------------------
print('Upsampling montly data...')