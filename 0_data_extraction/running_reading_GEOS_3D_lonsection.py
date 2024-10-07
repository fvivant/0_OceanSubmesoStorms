# Coded by Hector Torres (NASA-JPL)
# Adapted by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program for extracting 3D GEOS5 data at a given longitude
# used for the meridional cross sections of Fig. 5
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
varnames=['V','W','H','DTHDTCN'] # V, W, H, Convective diabatic heating
collections = [ 'inst_01hr_3d_'+i+'_Mv' for i in varnames ]
# varnames=['CNVMFC'] # Convective mass flux
# collections = [ 'tavg_01hr_3d_'+i+'_Mv' for i in varnames ]
# Note that for tavg_01hr_3d, M1 and M2 must be 30
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
lonsection=152.5
loni = [lonsection,lonsection+0.1]
lati = [31,41.85]

## Vertical levels
levels = np.arange(0,51,1)

## Data directory
diro = '/GEOSgcm_output/'
expid='c1440_llc2160'
diros=[diro+collections[i] for i in range(len(collections))]
flists=[  [datetime.strftime(n,diros[i]+'/DYAMOND_'+expid+'.'+collections[i]+'.' + '%Y%m%d_%H%M' + 'z.nc4') for n in date ] for i in range(len(diros))  ]

## Folder where data are saved
DIRC = '/path/to/your/extracted_data/2Dv/' #2Dv for "2D vertical section"
# i.e. meridional cross section

# ----------------- Extraction
print("Extraction/Interpolation will start")

for n in range(len(flists)):
    for i in range(0,len(flists[n])):
        ds,lon,lat,z = reading.GEOS_xr(flists[n][i],varnames[n],loni,lati,levels)
        name=varnames[n]+'_2Dv_'+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc'
        set=DATANC(varnames[n],lon,lat,z,ds)
        setf=set.sel(lon=lonsection)
        setf.to_netcdf(DIRC+'/'+name)
    print(varnames[n]+'  successfully extracted --> time spent = ',round((time.time()-t0)/60,2),' min')

print("SUCCESS if no error prompt --> total time : ",round((time.time()-t0)/60,2),' min',\
 " / time per variable : ", round((time.time()-t0)/60/len(varnames),2),' min')


