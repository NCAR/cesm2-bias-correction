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

from numba import vectorize, float64, jit, njit, guvectorize

import logging
logger = logging.getLogger(__name__)

# from distributed import Client
# from ncar_jobqueue import NCARCluster
# from pop_tools.grid import _compute_corners

def regrid_sst(in_tos, in_ta):
   """ Interpolate SST to CESM domain """

   logger.info("Interpolating SST to CESM domain")
   sst = in_tos.cf['surface_temperature']

   # Create a mask (not needed for interpolating to atmospheric grid, but just in case there are missing values)
   in_tos["mask"] = ~sst.cf.isel(time=0).isnull()

   # Regrid SST grid to whatever the atmospheric grid is automatically
   regrid = xe.Regridder(in_tos, in_ta, method = 'bilinear', periodic=True, unmapped_to_nan=True)
   regrid.to_netcdf('weights_gx1v6_latlon.nc') # write out weights for reuse 

   regridded_sst = regrid(in_tos)

   if logger.level == logging.DEBUG:
      print(regridded_sst)

   return regridded_sst 

def regrid_seaice(in_sic, in_ta):
   """ Interpolate SEAICE to CESM domain """

   logger.info("Interpolating SEAICE to CESM domain ...")
   ice_daily = in_sic['aice_d']*0.01
   in_sic["mask"] = ~ice_daily.isel(time=0).isnull()
   regrid = xe.Regridder(in_sic, in_ta, method = 'bilinear', periodic=True, unmapped_to_nan=True)

   regridded_ICE_DAILY = regrid(in_sic)

   return regridded_ICE_DAILY

@vectorize([float64(float64,float64,float64,float64)],nopython=True)
def pslec_atomic(temp_bot,phi_s,ps,pressure_bot): 
   """ Calculate sea level pressure (atomic execution using Numba vectorize) """

   LAPSE_RATE = 0.0065     #Kelvin per meter
   GRAV_CONST = 9.80616    #Meters per second per second
   SPEC_GAS_CONST = 287.04 #Joules per kilogram per Kelvin
   ALPHA_0 = LAPSE_RATE*SPEC_GAS_CONST/GRAV_CONST
    
   temp_surf = temp_bot*(1 + ALPHA_0*(ps/pressure_bot - 1)) #3b.5
   temp_bot_lapse = temp_surf + LAPSE_RATE*phi_s/GRAV_CONST #denoted T_0 in doc, 3b.13

    #These cases are partitions - there is no overlap in cases here. 
   if abs(phi_s/GRAV_CONST) < 1e-4: 
      psl = ps
   elif temp_surf >= 255 and temp_bot_lapse <= 290.5: 
      combo_term = ALPHA_0*phi_s/SPEC_GAS_CONST/temp_surf 
      psl =  ps*np.exp(combo_term/ALPHA_0*(1-1/2*(combo_term)+1/3*(combo_term)**2))
   elif temp_surf > 290.5 and temp_bot_lapse > 290.5: 
      T_star_modified = 1/2*(290.5+temp_surf) 
      psl = ps*np.exp(phi_s/SPEC_GAS_CONST/T_star_modified)
   elif temp_surf >=255 and temp_surf <= 290.5 and temp_bot_lapse > 290.5: 
      combo_term = 290.5-temp_surf
      psl = ps*np.exp(phi_s/SPEC_GAS_CONST/temp_surf*(1-1/2*(combo_term)+1/3*(combo_term)**2))
   elif temp_surf < 255 and temp_bot_lapse <= 290.5: 
      T_star_modified = 1/2*(255+temp_surf)
      combo_term = ALPHA_0*phi_s/SPEC_GAS_CONST/T_star_modified 
      psl =  ps*np.exp(combo_term/ALPHA_0*(1-1/2*(combo_term)+1/3*(combo_term)**2))
   elif temp_surf < 255 and temp_bot_lapse > 290.5: 
      alpha = SPEC_GAS_CONST/phi_s*(290.5-temp_surf)
      T_star_modified = 1/2*(255+temp_surf)
      combo_term = alpha*phi_s/SPEC_GAS_CONST/T_star_modified 
      psl = ps*np.exp(combo_term/alpha*(1-1/2*(combo_term)+1/3*(combo_term)**2))    
    
   return psl 
    
def pslec(temp_bottom: xr.DataArray, phi_surf: xr.DataArray, pressure_surf: xr.DataArray, pressure_bot: xr.DataArray): 
    return xr.apply_ufunc(pslec_atomic,temp_bottom,phi_surf,pressure_surf,pressure_bot)

def pres_on_hybrid():
   """ Calculate pressure on hybrid sigma-pressure levels """

