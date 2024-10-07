# Coded by Hector Torres (NASA-JPL)
# Adapted by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program required for extracting 3D GEOS5 data
# Requirements: 
## 1) llcmap_nea_split.py 
## 2) "geos_c1440_lats_lons_2D.nc" file available at:
## https://portal.nccs.nasa.gov/datashare/G5NR/DYAMONDv2/GEOS_6km_Atmosphere-MITgcm_4km_Ocean-Coupled/
## in the folder GEOSgcm_output/grid_coordinates
## 3) One must need to modify grd_dir by your own path

import numpy as np
import xarray as xr
from llcmap_nea_split import LLCMap_nea_split

def GEOS_xr(name,VAR,lon,lat,levels):
    dx = 0.07
    lat_out=np.arange(lat[0],lat[1]+dx,dx)
    lon_out=np.arange(lon[0],lon[1],dx)
  
    output=xr.DataArray(np.zeros((lat_out.shape[0],lon_out.shape[0])), \
                      coords=[ lat_out,lon_out], \
                      dims=[ 'lat', 'lon'])
    
    # grd_dir='/nobackup/fvivant/files/geos_c1440_lats_lons_2D.nc'
    grd_dir = "/GEOSgcm_output/grid_coordinates/geos_c1440_lats_lons_2D.nc"

    coords = xr.open_dataset(grd_dir)
    lon = coords.lons
    nx,ny = lon.shape

    XC=xr.where(coords.lons>=180,coords.lons-360,coords.lons)
    YC=coords.lats
    mapper = LLCMap_nea_split(YC.values,XC.values,lat_out,lon_out,radius=15e3)

    ## reading file
    ds1 = xr.open_mfdataset(name)
    
    ## empty matrix
    MAT = np.empty((len(levels),len(lat_out),len(lon_out)))
    z = np.empty((len(levels)))
    
    I = 0
    for k in levels:
        TMP=ds1[VAR][0,k,:,:].values
        output[:]=mapper((TMP))
        MAT[I,:,:] = output
        z[I] = ds1['lev'][k].values
        I += 1

    return MAT,lon_out,lat_out,z
