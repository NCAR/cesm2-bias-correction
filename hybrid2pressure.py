#!/usr/bin/env python3

""" 
This code mirrors the NCL script 'convert_cesm_hybric_nc_to_pressure_int.ncl', available at 
https://rda.ucar.edu/datasets/ds316.1/

The old NCL code bilinearly interpolates the grid, converts to pressure level coordinates, and converts to WRF-INT format.
This code will only interpolate the grid and convert to pressure level coordinates (no conversion to WRF-INT).

INPUTS:    3D CESM 6-hourly data, in netCDF4 format
OUTPUTS:   3D CESM 6-hourly data, interpolated onto lat/lon grid and in pressure coordinates, in netCDF4 format. 

Author: Kenton Wu 
Date: Jun 26, 2023 
Email: wukenton@ucar.edu, wukenton@gmail.com, wukenton@utexas.edu 
"""

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
import logging

# from distributed import Client
# from ncar_jobqueue import NCARCluster
# from pop_tools.grid import _compute_corners

def main(args):

   """ Read input data """
   logger.info("Opening files...")

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
   in_sic = xr.open_dataset("atmos_sic_1.nc")     # daily SEAICE percentage on POP grid (gaussian)

   """ Interpolate SST and SEAICE to CESM domain """
   logger.info('Converting Parallel Ocean Program data to coordinate system of atmospheric grid...')

   #Create a mask (not needed for interpolating to atmospheric grid, but just in case there are missing values)
   in_tos["mask"] = xr.where(~np.isnan(in_tos["tos"].sel(time=in_tos["tos"].time[0]),1,0)) 
   #in_ta["mask"] = xr.where(~np.isnan(in_ta["tos"].sel(time=in_ta["tos"].time[0]),1,0)) 

   #Regrids SST grid to whatever the atmospheric grid is automatically
   regrid = xe.Regridder(in_tos, in_ta, method = 'bilinear', periodic=True)
   regrid.to_netcdf('weights_gx1v6_latlon.nc') #write out weights for reuse 

   regridded_SST = regrid(in_tos)
   regridded_SST.to_netcdf('python_regrid_ncl.nc')

   #use some sort of broadcasting or view here to clone to a 6-hrly variable


   """ Interpolate monthly data to 6-hourly """
   print('Upsampling montly data...')

    

def regrid_sst(in_tos):
   """ Interpolate SST to CESM domain """

   logger.info("Interpolating SST to CESM domain")
   SST = in_tos.cf['surface_temperature']

   """ Create a mask (not needed for interpolating to atmospheric grid, but just in case there are missing values) """
   in_tos["mask"] = ~SST.cf.isel(time=0).isnull()

   # Regrid SST grid to whatever the atmospheric grid is automatically
   regrid = xe.Regridder(in_tos, in_ta, method = 'bilinear', periodic=True, unmapped_to_nan=True)
   regrid.to_netcdf('weights_gx1v6_latlon.nc') # write out weights for reuse 

   regridded_SST = regrid(in_tos)

   if logger.level == logging.DEBUG:
      print(regridded_SST)

   return regridded_SST   

def regrid_seaice(in_sic):
   """ Interpolate SEAICE to CESM domain """
   logger.info("Interpolating SEAICE to CESM domain ...")
   ICE_DAILY = in_sic['aice_d']*0.01
   in_sic["mask"] = ~ICE_DAILY.isel(time=0).isnull()
   regrid = xe.Regridder(in_sic, in_ta, method = 'bilinear', periodic=True, unmapped_to_nan=True)

   regridded_ICE_DAILY = regrid(in_sic)

   return regridded_ICE_DAILY

def configure_log(**kwargs):
   """ Configure logging """
   logpath = '/gpfs/csfs1/collections/rda/work/tcram/siparcs-2023/wukenton/cesm-bias-correction/logs'
   file = os.path.basename(__file__)
   logfile = '{}/{}.log'.format(logpath, os.path.splitext(file)[0])

   if 'loglevel' in kwargs:
      loglevel = kwargs['loglevel']
   else:
      loglevel = 'info'

   level = getattr(logging, loglevel.upper())
   format = '%(asctime)s - %(name)s - %(lineno)d - %(levelname)s - %(message)s'
   logging.basicConfig(filename=logfile, level=level, format=format)

   return

def parse_opts():
   """ Parse command line arguments """
   import argparse
   import textwrap
	
   desc = "Convert CESM data from hybrid sigma-pressure vertical coordinates to pressure level coordinates."	
   epilog = textwrap.dedent('''\
   Example:
        hybrid2pressure.py --case 20THC --output output_file.nc 
   ''')

   parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=desc, epilog=textwrap.dedent(epilog))
   parser.add_argument('-c', '--case', action="store", required=True, choices=['20THC', 'SSP126', 'SSP245', 'SSP370', 'SSP585'], help='CMIP6 historical or future scenario (20THC, SSP126, SSP245, SSP370, or SSP585)')
   parser.add_argument('-o', '--output', type=str, help='File name of desired output pressure levels')
   parser.add_argument('-l', '--loglevel', default="info", choices=['debug', 'info', 'warning', 'error', 'critical'], help='Set the logging level.  Default = info.')
   parser.add_argument('-w', '--weights', type=str, help="File name if reusing regridding weights")
   
   if len(sys.argv)==1:
      parser.print_help()
      sys.exit(1)

   args = parser.parse_args(sys.argv[1:])

   return args

if __name__ == "__main__":
   args = parse_opts()
   configure_log(loglevel=args.loglevel)
   logger = logging.getLogger(__name__)
   main(args)
