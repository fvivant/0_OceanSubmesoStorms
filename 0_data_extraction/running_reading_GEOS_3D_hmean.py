# Coded by Hector Torres (NASA-JPL)
# Adapted by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program for extracting 3D GEOS5 data averaged between two given heights
# used in Fig. 5c and 5f ( w(z<2km) )
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
varnames=['W'] # Convective mass flux and convective diabatic heating
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

## Vertical levels
levels = np.arange(0,51,1)

## Data directory
diro = '/GEOSgcm_output/'
expid='c1440_llc2160'
diros=[diro+collections[i] for i in range(len(collections))]
flists=[  [datetime.strftime(n,diros[i]+'/DYAMOND_'+expid+'.'+collections[i]+'.' + '%Y%m%d_%H%M' + 'z.nc4') for n in date ] for i in range(len(diros))  ]
# Height is needed for the interpolation:
diros_H=diro+'inst_01hr_3d_H_Mv'
flists_h=[datetime.strftime(n,diros_H+'/DYAMOND_'+expid+'.'+'inst_01hr_3d_H_Mv'+'.' + '%Y%m%d_%H%M' + 'z.nc4') for n in date ]

## Folder where data are saved
DIRC = '/path/to/your/extracted_data/2D/'

## Average boundary heights in meters
zlim1=0
zlim2=2000

# ----------------- Extraction
print("Extraction/Interpolation will start")

for n in range(len(flists)):
    for i in range(0,len(flists[n])):
        ds,lon,lat,z = reading.GEOS_xr(flists[n][i],varnames[n],loni,lati,levels)
        dh,lon,lat,z = reading.GEOS_xr(flists_h[n],'H',loni,lati,levels)
        set=DATANC(varnames[n],lon,lat,z,ds)
        setH=DATANC('H',lon,lat,z,dh)
        setf=set.where(setH.H>=zlim1).where(setH.H<=zlim2).mean(dim='depth',skipna=True)
        name=varnames[n]+'_'+str(zlim1)+'mean'+str(zlim2)+'_2D_'+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc'
        setf.to_netcdf(DIRC+'/'+name)

    print(varnames[n]+'  successfully extracted --> time spent = ',round((time.time()-t0)/60,2),' min')


print("SUCCESS if no error prompt --> total time : ",round((time.time()-t0)/60,2),' min',\
 " / time per variable : ", round((time.time()-t0)/60/len(varnames),2),' min')
