# Coded by Hector Torres (NASA-JPL)
# Adapted by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Requirements: 
## 1) reading_MITgcm_3D.py, llcmap_nea_split.py and face_connections.py 
## 2) data are available at:
## https://portal.nccs.nasa.gov/datashare/G5NR/DYAMONDv2/GEOS_6km_Atmosphere-MITgcm_4km_Ocean-Coupled/
## 3) xgcm must be version 0.6.1 !!! Issues with 0.8.1

import numpy as np
import xarray as xr
import reading_MITgcm_3D as reading

# ------------- Function
def timeline(initial_time,end_time,reference_time,
                 timdelta,rate,dt):

     """ Return timesteps corresponding to 
            the period selected
        Inputs:
        initial_time = '20200301000000' # YYYYMMDDHHMMSS
        end_time = '20200301000000'# YYYYMMDDHHMMSS
        reference_time (information in data.cal) = '20200119210000'
        timedelta = 1 hour
        dt = 45
        rate = 3600 # sampling rate (in seconds)
     """
     from datetime import datetime,timedelta
     from dateutil.parser import parse
     #from datetime import timedelta
     tim0 = parse(initial_time)
     tim1 = parse(end_time)
     tref = parse(reference_time)
     date=np.arange(tim0,tim1+timedelta(hours=timdelta),
                            timedelta(hours=timdelta)).astype(datetime)
     steps = int((tim0-tref).total_seconds()/dt)##self.rate
     i0 = steps
     del steps
     steps = int((tim1-tim0).total_seconds()/dt) ##self.rate
     i1 = i0+steps
     timesteps = np.arange(i0,i1+(rate/dt),rate/dt)
     return timesteps,date

# ----------------- Initialization

## Variable to extract
VAR = 'Theta' # 'U','V','W','Eta'

## Dates
# initial_time = '20201201000000'
# end_time = '20210228000000'
# test for one date:
initial_time = '20201201000000'
end_time = '20201201000000'

#### don't touch ############
reference_time = '20200119210000'
Timedelta  = 1 #hour
rate = 3600 # seconds
dt = 45 # seconds 
#############################
all_iters,date=timeline(initial_time,end_time,reference_time,
                 Timedelta,rate,dt)
# ---------------------------------

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
levels = [0] # surface level of the MITgcm, one can extract more levels

## Folder where data are saved
diro = '/path/to/your/extracted_data/2D/'
# note that one need to modify the data directory to extract (diro) in reading_MITgcm_3D.py

for i in range(0,len(all_iters)):
     data,lon,lat,zaxis = reading.getting_3D(VAR,
                    all_iters[i],
                                loni,lati,levels)
     nam = VAR+'_3D_'+date[i].strftime("%Y%m%d_%H%M")+'_ocean.nc'

     # data2,lon2,lat2,zaxis2 = reading.getting_3D(VAR,
     #                all_iters[i],
     #                            loni2,lati2,levels)

     # Here, 2D surface values are stored (e.g, the SST)
     # but one could add the depth dimension bellow and store "data[:]"
     output=xr.Dataset(coords={'lat':lat,'lon':lon,})
     output[VAR] = (('lat','lon'),data[0]) 

     # output2=xr.DataArray(np.zeros((zaxis2.shape[0],lat2.shape[0],lon2.shape[0])),
     #                  dims=("depth","lat","lon"),
     #                  coords={'lat':lat2,'lon':360+lon2,
     #                         'depth':zaxis2})

     # output2[VAR] = (['depth','lat','lon'],data2)

     outputf=output #xr.concat([output,output2],dim="lon")
     # print("store_start")
     outputf.to_netcdf(diro+'/'+nam)
     # print("store_end")

