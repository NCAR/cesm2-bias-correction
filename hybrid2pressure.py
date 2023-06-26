#This code mirrors the ncl script convert_cesm_hybric_nc_to_pressure_int.ncl
#The old code interpolates the grid, convert to pressure coordinates, and converts to WRF-INT format
#This code will only interpolate the grid and convert to pressure coordinates (no conversion to WRF-INT)

#INPUTS:    3D CESM 6hr data, in netCDF4 format
#OUTPUTS:   3D CESM 6hr data, interpolated onto lat/lon grid and in pressure coordinates, in netCDF4 format. 

#Author: Kenton Wu 
#Date: Jun 26, 2023 
#Email: wukenton@ucar.edu, wukenton@gmail.com, wukenton@utexas.edu 

import argparse 
import xarray as xr 


#Command line option handling ----------------------------------------------------------------------------------
parser = argparse.ArgumentParser()  

parser.add_argument('CASE',type=str, help='One of the following IPCC Climate Scenarios: 20THC/RCP85/RCP60/RCP45')
parser.add_argument('--o',type=str,help='Output directory path')

#File Handling ----------------------------------------------------------------------------------

print("Opening files...")
in_ta = open("atmos_ta.nc", "r")         # 6-hourly 3-d T
in_ua = open("atmos_ua.nc", "r")         # 6-hourly 3-d U
in_va = open("atmos_va.nc", "r")         # 6-hourly 3-d V
in_hus = open("atmos_hus.nc", "r")       # 6-hourly 3-d Q
in_ps = open("atmos_ps.nc", "r")         # 6-hourly surface pressure
in_zsfc = open("atmos_zsfc.nc", "r")     # static surface geopotential
in_lmask = open("atmos_lmask.nc", "r")   # static land mask
in_snw = open("atmos_snw_1.nc", "r")     # monthly SWE
in_mrlsl = open("atmos_mrlsl_1.nc", "r") # monthly soil moisture
in_ts = open("atmos_ts_1.nc", "r")       # monthly skin temp
in_tsl = open("atmos_tsl_1.nc", "r")     # monthly soil temp
in_tos = open("atmos_tos_1.nc", "r")     # daily SST on pop grid (gaussian)
in_sic = open("atmos_sic_1.nc", "r")     # daily SEAICE % on POP grid (gaussian)

#Read variables into file 