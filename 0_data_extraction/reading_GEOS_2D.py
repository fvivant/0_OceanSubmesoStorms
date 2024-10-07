# Coded by Hector Torres (NASA-JPL)
# Adapted by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program required for extracting 2D GEOS5 data
# Requirements: 
## 1) llcmap_nea_split.py 
## 2) "geos_c1440_lats_lons_2D.nc" file available at:
## https://portal.nccs.nasa.gov/datashare/G5NR/DYAMONDv2/GEOS_6km_Atmosphere-MITgcm_4km_Ocean-Coupled/
## in the folder GEOSgcm_output/grid_coordinates
## 3) One must need to modify GEOS_gridfile by your own path

import numpy as np
import xarray as xr
from llcmap_nea_split import LLCMap_nea_split


def GEOS_xr_coll_date_location_fol(filename, VAR, lat1, lat2, lon1, lon2):
  dx=0.07
  lat_out = np.arange(lat1, lat2 + dx, dx)
  lon_out = np.arange(lon1, lon2 , dx) 

  output=xr.DataArray(np.zeros((lat_out.shape[0],lon_out.shape[0])), \
                      coords=[ lat_out,lon_out], \
                      dims=[ 'lat', 'lon'])


  GEOS_gridfile = "/GEOSgcm_output/grid_coordinates/geos_c1440_lats_lons_2D.nc"
  
  gridds = xr.open_dataset(GEOS_gridfile)
  XC = gridds.lons
  XC = xr.where(XC>=180, XC- 360, XC)
  YC = gridds.lats

  mapper = LLCMap_nea_split(YC.values, XC.values,lat_out,lon_out,radius=15e3)

  ds1 = xr.open_dataset(filename)
    
  TMP=ds1[VAR]
    
  output[:] = mapper(TMP.mean('time').values)
  map = output
  return lon_out,lat_out,map

    
