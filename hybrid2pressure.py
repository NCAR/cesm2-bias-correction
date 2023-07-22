#This code mirrors the ncl script convert_cesm_hybric_nc_to_pressure_int.ncl
#The old code bilinearly interpolates the grid, convert to pressure coordinates, and converts to WRF-INT format
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

print('Converting Parallel Ocean Program data to coordinate system of atmospheric grid...')

#Create a mask (not needed for interpolating to atmospheric grid, but just in case there are missing values)
in_tos["mask"] = xr.where(~np.isnan(in_tos["tos"].sel(time=in_tos["tos"].time[0]),1,0)) 
#in_ta["mask"] = xr.where(~np.isnan(in_ta["tos"].sel(time=in_ta["tos"].time[0]),1,0)) 

#Regrids SST grid to whatever the atmospheric grid is automatically
regrid = xe.Regridder(in_tos, in_ta, method = 'bilinear', periodic=True)
regrid.to_netcdf('weights_gx1v6_latlon.nc') #write out weights for reuse 

regridded_SST = regrid(in_tos)
regridded_SST.to_netcdf('python_regrid_ncl.nc')
#
#use some sort of broadcasting or view here to clone to a 6-hrly variable


#Turn monthly data into 6hr data ----------------------------------------------------------------------------------
print('Upsampling montly data...')

def is_missing(data : xr.DataArray): 
    