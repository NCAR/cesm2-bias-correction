#!/usr/bin/env python3

"""
Script to process CESM data
1. Interpolate SST and SEAICE to CESM atmospheric domain
2. Interpolate monthly variables to 6-hourly
3. Interpolate 3-d variables to vertical pressure levels
4. Write NetCDF output
"""

import logging
import logging.handlers
import os, sys
import xarray as xr
from cdo import *
from pathlib import Path
from datetime import datetime

from cesm_bias_correction.hybrid2pressure import *

cases = {
   '20THC': {
      'start_year': 1850,
      'end_year': 2014
   },
   'RCP45': {
      'start_year': 2006,
      'end_year': 2100
   },
   'RCP60': {
      'start_year': 2006,
      'end_year': 2100
   },
   'RCP85': {
      'start_year': 2006,
      'end_year': 2100
   },
   'SSP126': {
      'start_year': 2015,
      'end_year': 2100
   },
   'SSP245': {
      'start_year': 2015,
      'end_year': 2100
   },
   'SSP370': {
      'start_year': 2015,
      'end_year': 2100
   },
   'SSP585': {
      'start_year': 2015,
      'end_year': 2100
   }
}

_DAYS_IN_MONTH = [31,28,31,30,31,30,31,31,30,31,30,31]  # Ignore leap days

def main(opts):

   # iterate through years and months
   for year in range(opts['start_year'], opts['end_year']):
      for month in range(1, 13):
         

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

   # sea level pressure
   slp = pslec(temp.isel(lev=-1), phi_surf, surf_pressure, P_hybrid.isel(lev=-1))

   # pressure on hybrid sigma-pressure levels
   P_hybrid = pres_on_hybrid_ccm(surf_pressure,hyam,hybm).astype(np.single)

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
   """ Parse command line arguments.  Returns a dictionary of command line arguments and values. """

   import argparse
   import textwrap
	
   desc = "Convert CESM data from hybrid sigma-pressure vertical coordinates to pressure level coordinates."	
   epilog = textwrap.dedent('''\
   Example:
        process_cesm.py --case 20THC --output output_file.nc 
   ''')

   parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=desc, epilog=textwrap.dedent(epilog))
   parser.add_argument('-c', '--case', action="store", required=True, choices=['20THC', 'SSP126', 'SSP245', 'SSP370', 'SSP585'], help='CMIP6 historical or future scenario (20THC, SSP126, SSP245, SSP370, or SSP585)')
   parser.add_argument('-s', '--start-year', action="store", required=True, type=int, help="Start year to process")
   parser.add_argument('-e', '--end-year', action="store", help="End year to process")
   parser.add_argument('-o', '--output', type=str, help='File name of desired output pressure levels')
   parser.add_argument('-w', '--weights', type=str, help="Input file name of regridding weights if reusing from a previous run.")
   parser.add_argument('-l', '--loglevel', default="info", choices=['debug', 'info', 'warning', 'error', 'critical'], help='Set the logging level.  Default = info.')
   
   if len(sys.argv)==1:
      parser.print_help()
      sys.exit(1)

   args = parser.parse_args(sys.argv[1:])
   opts = vars(args)

   case = opts['case']
   case_start_year = cases[case]['start_year']
   case_end_year = cases[case]['end_year']

   if opts['end_year'] is None:
      logger.info('Setting end year to same as start-year {}.'.format(opts['start_year']))
      opts['end_year'] = opts['start_year']

   if opts['start_year'] > case_end_year:
      print("Case {}: start year must be earlier than {}". format(case, case_end_year))
      sys.exit(1)
   if opts['start_year'] < case_start_year:
      print("Case {}: start year must be equal or greater than {}.  Setting to {}.".format(case, case_start_year, case_start_year))
      opts['start_year'] = case_start_year
   if opts['end_year'] < case_start_year or opts['end_year'] < opts['start_year']:
      print("Case {}: end year must be equal or greater than {}.  Setting to {}.".format(case, case_start_year, case_start_year))
      opts['end_year'] = opts['start_year']

   logger.info('Processing CESM for case {}, years {} to {}'.format(opts['case'], opts['start_year'], opts['end_year']))

   return opts

if __name__ == "__main__":
   opts = parse_opts()
   configure_log(loglevel=opts['loglevel'])
   logger = logging.getLogger(__name__)
   main(opts)
