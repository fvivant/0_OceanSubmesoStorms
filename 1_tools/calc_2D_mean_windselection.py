# Coded by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program for computing the mean of 2D fields with northward or southward wind section
# used in Fig. 4, Fig. 6, Extended Fig. 4, Extended Fig. 7

import xarray as xr
import numpy as np
from datetime import datetime,timedelta
import time

# ---------------- Function
def Mean_wind(dates,name,varname,folder,folderout,V10M_mean,mode,lim,date_name):
        ## Time mean with wind section (all or northwards or southwards)
        
        # Handle different name-formats of variables 
        # In particular, MITgcm variables have been stored with '_ocean' at the end
        if name=='Thetaoc':
                name='Theta'
                file=datetime.strftime(dates[0],folder+name+'_3D_'+'%Y%m%d_%H%M' + '_ocean.nc')
                var0=xr.open_dataset(file)
                lonn,latt=var0.lon.data,var0.lat.data
                llat,llon=var0[varname].shape
                lt=len(dates)
                
                var=xr.Dataset( {varname: ( ('time','lat','lon'), np.zeros((lt,llat,llon)) )},\
                        coords={'time': dates,"lat": latt, "lon": lonn } )
                
                for i in range(len(dates)):
                        file=datetime.strftime(dates[i],folder+name+'_3D_'+'%Y%m%d_%H%M' + '_ocean.nc')
                        var[varname][i,:]=xr.open_dataset(file)[varname].data
        
                        
        else:
                file=datetime.strftime(dates[0],folder+name+'_2D_'+'%Y%m%d_%H%M' + '.nc')
                var0=xr.open_dataset(file)
                lonn,latt=var0.lon.data,var0.lat.data
                llat,llon=var0[varname].shape
                lt=len(dates)
                
                var=xr.Dataset( {varname: ( ('time','lat','lon'), np.zeros((lt,llat,llon)) )},\
                        coords={'time': dates,"lat": latt, "lon": lonn } )
                
                for i in range(len(dates)):
                        file=datetime.strftime(dates[i],folder+name+'_2D_'+'%Y%m%d_%H%M' + '.nc')
                        var[varname][i,:]=xr.open_dataset(file)[varname].data
        
        ## NW for northwards, SW for southwards, all for every winds
        if mode=='NW':
                var[varname].where(V10M_mean>=0)\
                        .mean(dim='time').to_netcdf(folderout+name+'_mean_'+date_name+'_NW'+str(abs(lim))+'sup_2D.nc')#
        elif mode=='SW':
                var[varname].where(V10M_mean<=lim)\
                        .mean(dim='time').to_netcdf(folderout+name+'_mean_'+date_name+'_SW'+str(abs(lim))+'_2D.nc')
        elif mode=='all':
                var[varname].mean(dim='time').to_netcdf(folderout+name+'_mean_'+date_name+'_2D.nc')
        
        
# ---------------- Execution

### 1 - Mean meridonal wind calculation

# 5-day mean (Fig. 4 and 6)
datefull = np.arange(datetime(2020,12,16,0),datetime(2020,12,21,0),timedelta(hours=1))\
       .astype(datetime)  
date_name='16dec21dec' # storage name

# 2-month mean
# datefull = np.arange(datetime(2020,12,2,0),datetime(2021,1,31,0),timedelta(hours=1))\
#        .astype(datetime)
# date_name='decjan' # storage name
 
dates=datefull
 
# Input data folder:      
folder='/path/to/your/extracted_data/2D/'

name,varname='V10M','V10M'
# getting shape
file=datetime.strftime(dates[0],folder+name+'_2D_'+'%Y%m%d_%H%M' + '.nc')
var0=xr.open_dataset(file)
lonn,latt=var0.lon.data,var0.lat.data
llat,llon=var0[varname].shape
lt=len(dates)
# creating dataset
var=xr.Dataset( {varname: ( ('time','lat','lon'), np.zeros((lt,llat,llon)) )},\
        coords={'time': dates,"lat": latt, "lon": lonn } )

# getting merdional wind values
for i in range(len(dates)):
        file=datetime.strftime(dates[i],folder+name+'_2D_'+'%Y%m%d_%H%M' + '.nc')
        var[varname][i,:]=xr.open_dataset(file)[varname].data
        
# mean meridional wind
V10M_mean=var[varname].mean(dim='lon').mean(dim='lat')


### 2 - Time mean with wind selection calculation

# Field of Fig. 4
# Same for Fig. 6 but with 'DTHDTCN_2000m'
listname=['U10M','V10M','Thetaoc','EFLUX','CAPE','PBLH','gradTheta_dwind','divU10','MFC','DTHDTCN_3000m','PRECCON']
listvar=['U10M','V10M','Theta','EFLUX','CAPE','PBLH','gradTheta_dwind','divU10','MFC','DTHDTCN','PRECCON']

# Output folder:
folderout='/path/to/your/extracted_data/2D/'

for i in range(len(listname)):
        t0=time.time()
        Mean_wind(datefull,listname[i],listvar[i],folder,folderout,
                  V10M_mean,'NW',0,date_name) # NW  for northward winds
        print(time.time()-t0)
# -------------------------------      

## Special handling for CNVMFC, which is stored at hh_30
# Field of Fig. 4
# Same for Fig. 6 but with 'CNVMFC_2000m'
listname=['CNVMFC_3000m']
listvar=['CNVMFC']

# Output folder:
folderout='/path/to/your/extracted_data/2D/'

# Since V10M and CNVMFC are not stored at the same dates (30 min),
# we assume that we can perform the wind selection with a 30 min shift
shifted_dates=dates+timedelta(minutes=30)
V10M_mean=V10M_mean.assign_coords(time= shifted_dates)

for i in range(len(listname)):
        t0=time.time()
        Mean_wind(shifted_dates,listname[i],listvar[i],folder,folderout,
                  V10M_mean,'NW',0,date_name) # NW  for northward winds (SW  for southward winds)
        print(time.time()-t0)