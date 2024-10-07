# Coded by Hector Torres (NASA-JPL)
# Adapted by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program for extracting 2D GEOS5 data
# Requirements: 
## 1) reading_GEOS_2D.py and llcmap_nea_split.py 
## 2) data are available at:
## https://portal.nccs.nasa.gov/datashare/G5NR/DYAMONDv2/GEOS_6km_Atmosphere-MITgcm_4km_Ocean-Coupled/
## in the folder GEOSgcm_output/

import numpy as np
import reading_GEOS_2D as read
from datetime import datetime,timedelta
import xarray as xr
import time

# ------------- Functions
def timeline(coll,y1, m1, d1,h1, M1, y2, m2, d2,h2, M2):
    # Timeline maker for 2D GEOS data depending on the collection consider
    collection1='geosgcm_surf'
    # Hourly surface quantities collection. Contains for exemple T10M and Q10M.
    collection2 = 'const_2d_asm_Mx'
    # Constant repository - stores constants like FRLAND, FRLANDICE, FROCEAN, etc
    collection3 = 'inst_15mn_2d_asm_Mx'
    # Instanteneous 2D quantities collection stored every 15min
    # Contains in particular CAPE, U10M, V10M, QV2M (specific humidity at 2 meters)
    collection4 = 'tavg_15mn_2d_flx_Mx'
    # Averaged over 15min 2D quantities stored every 15min at hh_7min ,hh_22min, hh_37min and hh_52min
    # Contains in particular EFLUX (LHF), PBLH, PRECCON (convective precipitation), PRECTOT (total precipitation)
    
    if (coll == 'SURF'):
        collection = collection1
        M1 = 30
        M2 = 30
        deltatime = timedelta(hours=1)
    elif (coll == 'CONST'):
        collection = collection2
        y1 = 2020
        y2 = 2020
        m1 = 1
        m2 = 1
        d1 = 22
        d2 = 22
        M1 = 0
        M2 = 0
        h1 = 0
        h2 = 1
        deltatime = timedelta(hours=1)
    elif (coll == 'INST15'):
        collection = collection3
        deltatime = timedelta(hours=1)#timedelta(minutes=15)
        # Here, data are extracted every hours, one can choose every 15 minutes
    elif (coll == 'TAVG15MIN'):
        collection = collection4
        deltatime = timedelta(hours=1)#timedelta(minutes=15)
        # Here, data are extracted every hours, one can choose every 15 minutes
        M1 = 7 
        M2 = 7

    date = np.arange(datetime(y1,m1,d1,h1,M1),
                 datetime(y2,m2,d2,h2,M2), deltatime).astype(datetime)
    return date,collection


def DATANC2D(VAR,LON,LAT,DATA,poslon=False):
    #poslon==True when LON<0 for LON continuity when concat
    if poslon==True:
        output=xr.Dataset(coords={'lat':LAT,'lon':LON+360})

    else:
        output=xr.Dataset(coords={'lat':LAT,'lon':LON})

    output[VAR] = (('lat','lon'),DATA)

    return output

# ----------------- Initialization
## Variables to extract
# Two collections has been used in the manuscript
# coll = 'TAVG15MIN' # Contains in particular EFLUX (LHF), PRECCON, PRECTOT, PBLH
# varnames=['EFLUX','PBLH','PRECCON','PRECTOT']

coll = 'INST15'  # Contains in particular CAPE, U10M, V10M, QV2M
varnames=['CAPE','U10M','V10M','QV2M']
# ------------------

## Dates
# y1,m1,d1,h1,M1 = 2020,12,1,0,0 #start
# y2,m2,d2,h2,M2 = 2021,2,28,0,0 #end
# test for one date:
y1,m1,d1,h1,M1 = 2020,12,1,0,0 #start
y2,m2,d2,h2,M2 = 2020,12,1,1,0 #end
date,collection = timeline(coll,y1, m1, d1,h1, M1, y2, m2, d2,h2, M2)
# -------------------

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

## Data directory
diro = '/GEOSgcm_output/'
expid='c1440_llc2160'

## Dates in GEOS5 format
flist=[datetime.strftime(n,diro+'/DYAMOND_'+expid+'.'+collection+'.' + '%Y%m%d_%H%M' + 'z.nc4')\
        for n in date ]

## Folder where data are saved
DIRC = '/path/to/your/extracted_data/2D/'

# ----------------- Extraction
print("Extraction/Interpolation will start")

for i in range(len(varnames)):
    t0=time.time()
    for n in range(len(flist)):
        lon,lat,ds = read.GEOS_xr_coll_date_location_fol(flist[n],
                                                varnames[i],lati[0],
                                                lati[1], loni[0],
                                                loni[1])

        # lon2,lat2,ds2 = read.GEOS_xr_coll_date_location_fol(coll,flist[n],
        #                                         varnames[i],lati2[0],
        #                                         lati2[1], loni2[0],
        #                                         loni2[1])

        if coll=='TAVG15MIN':
            # -timedelta(minutes=7) to store data at hh_00, hh+1_00, hh+2_00, etc
            name=varnames[i]+'_2D_'+datetime.strftime(date[n]-timedelta(minutes=7),'%Y%m%d_%H%M')+'.nc'
        else:
            name=varnames[i]+'_2D_'+datetime.strftime(date[n],'%Y%m%d_%H%M')+'.nc'
            
        set=DATANC2D(varnames[i],lon,lat,ds.data)
        # set2=DATANC2D(varnames[i],lon2,lat2,ds2.data,poslon=True)
        setf=set #xr.concat([set,set2],dim="lon")
        setf.to_netcdf(DIRC+'/'+name)

    print(varnames[i]+'  successfully extracted --> time spent = ',round((time.time()-t0)/60,2),' min')


print("SUCCESS if no error prompt --> total time : ",round((time.time()-t0)/60,2),' min',\
 " / time per variable : ", round((time.time()-t0)/60/len(varnames),2),' min')

