# Coded by Hector Torres (NASA-JPL)
# Adapted by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program for extracting full 3D GEOS5 data
# The results presented in the manuscript does not require this program.
# Instead, running_reading_GEOS_3D_hinterp.py, running_reading_GEOS_3D_hmean.py
# and running_reading_GEOS_3D_lonsection.py have been used for the sake of storage.
# Those three programs are based on this one "running_reading_GEOS_3D.py"
# Requirements: 
## 1) reading_GEOS_3D.py and llcmap_nea_split.py 
## 2) data are available at:
## https://portal.nccs.nasa.gov/datashare/G5NR/DYAMONDv2/GEOS_6km_Atmosphere-MITgcm_4km_Ocean-Coupled/
## in the folder GEOSgcm_output/

from datetime import datetime,timedelta
import numpy as np
import xarray as xr
import reading_GEOS_3D as reading
import time

# ------------- Function
def DATANC(VAR,LON,LAT,Z,DATA,poslon=False):
    #poslon==True when LON<0 for LON continuity when concatanation
    if poslon==True:
        output=xr.Dataset(coords={'lat':LAT,'lon':LON+360,
                              'depth':Z})

    else:
        output=xr.Dataset(coords={'lat':LAT,'lon':LON,
                              'depth':Z})

    output[VAR] = (('depth','lat','lon'),DATA)

    return output

# ----------------- Initialization
t0=time.time()

## Variables to extract
varnames=['P'] #'U','V','W'
collections = [ 'inst_01hr_3d_'+i+'_Mv' for i in varnames ]
# ------------------

## Dates
# y1,m1,d1,h1,M1 = 2020,12,1,0,0 #start
# y2,m2,d2,h2,M2 = 2021,2,28,0,0 #end
# test for one date:
y1,m1,d1,h1,M1 = 2020,12,1,0,0 #start
y2,m2,d2,h2,M2 = 2020,12,1,1,0 #end
date = np.arange(datetime(y1,m1,d1,h1,M1),datetime(y2,m2,d2,h2,M2),timedelta(hours=1)).astype(datetime)
t = date.shape
# ------------------

## Domain
# Study domain:
loni=[142.3,162]
lati=[31,41.85]
# Domain in Figure 2a:
# loni=[130,170]
# lati=[25,50]
# ---------------------------------
# One can extract a larger domain with positive and negative longitudes
# To do so, one need to extract two domains because of the change
# of longitude sign in the North Pacific, and then concat them
# One can uncomment every lines bellow with loni2, lati2 to do so

# loni=[130,180]
# lati=[25,50]
# loni2 = [-180,-160]
# lati2=lati
# ---------------------------------

## Vertical levels
levels = np.arange(0,51,1)

## Data directory
diro = '/GEOSgcm_output/'
expid='c1440_llc2160'
diros=[diro+collections[i] for i in range(len(collections))]
flists=[  [datetime.strftime(n,diros[i]+'/DYAMOND_'+expid+'.'+collections[i]+'.' + '%Y%m%d_%H%M' + 'z.nc4') for n in date ] for i in range(len(diros))  ]

## Folder where data are saved
DIRC = '/path/to/your/extracted_data/3D/'

# ----------------- Extraction
print("Extraction/Interpolation will start")

for n in range(len(flists)):
    for i in range(0,len(flists[n])):
        ds,lon,lat,z = reading.GEOS_xr(flists[n][i],varnames[n],loni,lati,levels)
        # ds2,lon2,lat2,z2 = reading.GEOS_xr(flists[n][i],varnames[n],loni2,lati2,levels) 
        name=varnames[n]+'_3D_'+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc'
        set=DATANC(varnames[n],lon,lat,z,ds)
        # set2=DATANC(varnames[n],lon2,lat2,z2,ds2,poslon=True)
        setf=set #xr.concat([set,set2],dim="lon")
        setf.to_netcdf(DIRC+'/'+name)
    print(varnames[n]+'  successfully extracted --> time spent = ',round((time.time()-t0)/60,2),' min')

print("SUCCESS if no error prompt --> total time : ",round((time.time()-t0)/60,2),' min',\
 " / time per variable : ", round((time.time()-t0)/60/len(varnames),2),' min')

