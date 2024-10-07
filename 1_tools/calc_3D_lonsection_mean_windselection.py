# Coded by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program for computing the mean of a meridional cross section with northward or southward wind section
# used in Fig. 5, Extended Fig. 5 and Extended Fig. 6

import xarray as xr
import numpy as np
from datetime import datetime,timedelta


def mean(dates,folder,folder2D,folderout,lonsel,condi,condi_CNVMFC,savename):

        latmin,latmax=31,41.85
        
        # getting shape and coordinates
        name='V'
        path=folder+name+'_2Dv_'
        file=[path+datetime.strftime(date,'%Y%m%d_%H%M' + '.nc') for date in dates]
        V0=xr.open_dataset(file[0]).sel(lat=slice(latmin,latmax))
        ldepth,llat=V0.V.shape
        lt=len(dates)
        
        ############### 1 - # 3D fields for Fig. 5a, d
                                 
        # Height
        name='H'
        path=folder+name+'_2Dv_'
        file=[path+datetime.strftime(date,'%Y%m%d_%H%M' + '.nc') for date in dates]
        H=xr.Dataset( {'H': ( ('time','depth','lat'), np.zeros((lt,ldepth+1,llat)) )},\
                coords={'time': dates,"lat": V0.lat.data, "depth": np.array(list(V0.depth.data)+[72]) } )
        for i in range(len(dates)):
                H0=xr.open_dataset(file[i]).sel(lat=slice(latmin,latmax)).isel(lat=slice(0,llat)) 
                H.H[i,:-1,:]=H0.H.data
                H.H[i,-1,:]=np.ones(llat)*0 #2 meters surface layer
        path=folderout+name+'_2Dv_'
        H.H.where(condi).mean(dim='time',skipna=True).to_netcdf(path+savename)
        
        # Convective mass flux
        name='CNVMFC'
        CNVMFC=xr.Dataset( {'CNVMFC': ( ('time','depth','lat'), np.zeros((lt,ldepth,llat)) )},\
                coords={'time': dates+timedelta(minutes=30),"lat": V0.lat.data, "depth": np.array(list(V0.depth.data)) } )
        for i in range(len(dates)):
                file=datetime.strftime(dates[i]+timedelta(minutes=30),folder+name+'_2Dv_'+'%Y%m%d_%H%M' + '.nc')
                CNVMFC.CNVMFC[i,:]=xr.open_dataset(file).sel(lat=slice(latmin,latmax)) .CNVMFC.data
        path=folderout+name+'_2Dv_'
        CNVMFC.CNVMFC.where(condi_CNVMFC).mean(dim='time',skipna=True).to_netcdf(path+savename)
        
        # Convective diabatic heating
        name='DTHDTCN'
        DTHDTCN=xr.Dataset( {'DTHDTCN': ( ('time','depth','lat'), np.zeros((lt,ldepth,llat)) )},\
                coords={'time': dates,"lat": V0.lat.data, "depth": np.array(list(V0.depth.data)) } )
        for i in range(len(dates)):
                file=datetime.strftime(dates[i],folder+name+'_2Dv_'+'%Y%m%d_%H%M' + '.nc')
                DTHDTCN.DTHDTCN[i,:]=xr.open_dataset(file).sel(lat=slice(latmin,latmax)) .DTHDTCN.data
        path=folderout+name+'_2Dv_'
        DTHDTCN.DTHDTCN.where(condi).mean(dim='time',skipna=True).to_netcdf(path+savename)
        
        # Vertical velocity
        name='W'
        W=xr.Dataset( {'W': ( ('time','depth','lat'), np.zeros((lt,ldepth,llat)) )},\
                coords={'time': dates,"lat": V0.lat.data, "depth": np.array(list(V0.depth.data)) } )
        for i in range(len(dates)):
                file=datetime.strftime(dates[i],folder+name+'_2Dv_'+'%Y%m%d_%H%M' + '.nc')
                W.W[i,:]=xr.open_dataset(file).sel(lat=slice(latmin,latmax)) .W.data

        path=folderout+name+'_2Dv_'
        W.W.where(condi).mean(dim='time',skipna=True).to_netcdf(path+savename)
        
        # Meridonal velocity
        name='V'
        V=xr.Dataset( {'V': ( ('time','depth','lat'), np.zeros((lt,ldepth,llat)) )},\
                coords={'time': dates,"lat": V0.lat.data, "depth": np.array(list(V0.depth.data)) } )
        for i in range(len(dates)):
                file=datetime.strftime(dates[i],folder+name+'_2Dv_'+'%Y%m%d_%H%M' + '.nc')
                V.V[i,:]=xr.open_dataset(file).sel(lat=slice(latmin,latmax)).V.data
        path=folderout+name+'_2Dv_'
        V.V.where(condi).mean(dim='time',skipna=True).to_netcdf(path+savename)
        
        # PBLH (a 2D field)
        name='PBLH'
        PBLH=xr.Dataset( {'PBLH': ( ('time','lat'), np.zeros((lt,llat)) )},\
                coords={'time': dates,"lat": V0.lat.data } )
        for i in range(len(dates)):
                file=datetime.strftime(dates[i],folder2D+name+'_2D_'+'%Y%m%d_%H%M' + '.nc')
                PBLH.PBLH[i,:]=xr.open_dataset(file).sel(lon=lonsel,method='nearest').sel(lat=slice(latmin,latmax)).PBLH.data
        path=folderout+name+'_2Dv_'
        PBLH.PBLH.where(condi).mean(dim='time',skipna=True).to_netcdf(path+savename)
        
        ############### 2 - # 2D fields for Fig. 5b, c, e, f
        # SST
        name='Theta'
        SST=xr.Dataset( {'Theta': ( ('time','lat'), np.zeros((lt,llat)) )},\
                coords={'time': dates,"lat": V0.lat.data } )
        for i in range(len(dates)):
                file=datetime.strftime(dates[i],folder2D+name+'_3D_'+'%Y%m%d_%H%M' + '_ocean.nc')
                SST.Theta[i,:]=xr.open_dataset(file).sel(lon=lonsel,method='nearest').sel(lat=slice(latmin,latmax)).Theta.data
        path=folderout+name+'_ocean_2Dv_'
        SST.Theta.where(condi).mean(dim='time',skipna=True).to_netcdf(path+savename)
        
        # LHF                
        name='EFLUX'
        EFLUX=xr.Dataset( {'EFLUX': ( ('time','lat'), np.zeros((lt,llat)) )},\
                coords={'time': dates,"lat": V0.lat.data } )
        for i in range(len(dates)):
                file=datetime.strftime(dates[i],folder2D+name+'_2D_'+'%Y%m%d_%H%M' + '.nc')
                EFLUX.EFLUX[i,:]=xr.open_dataset(file).sel(lon=lonsel,method='nearest').sel(lat=slice(latmin,latmax)).EFLUX.data
                
        path=folderout+name+'_2Dv_'
        EFLUX.EFLUX.where(condi).mean(dim='time',skipna=True).to_netcdf(path+savename)
        
        # CAPE
        name='CAPE'
        CAPE=xr.Dataset( {'CAPE': ( ('time','lat'), np.zeros((lt,llat)) )},\
                coords={'time': dates,"lat": V0.lat.data } )
        for i in range(len(dates)):
                file=datetime.strftime(dates[i],folder2D+name+'_2D_'+'%Y%m%d_%H%M' + '.nc')
                CAPE.CAPE[i,:]=xr.open_dataset(file).sel(lon=lonsel,method='nearest').sel(lat=slice(latmin,latmax)).CAPE.data
        path=folderout+name+'_2Dv_'
        CAPE.CAPE.where(condi).mean(dim='time',skipna=True).to_netcdf(path+savename)
        
        # Convective precipitation
        name='PRECCON'
        PRECCON=xr.Dataset( {'PRECCON': ( ('time','lat'), np.zeros((lt,llat)) )},\
                coords={'time': dates,"lat": V0.lat.data } )
        for i in range(len(dates)):
                file=datetime.strftime(dates[i],folder2D+name+'_2D_'+'%Y%m%d_%H%M' + '.nc')
                PRECCON.PRECCON[i,:]=xr.open_dataset(file).sel(lon=lonsel,method='nearest').sel(lat=slice(latmin,latmax)).PRECCON.data
        path=folderout+name+'_2Dv_'
        PRECCON.where(condi).mean(dim='time',skipna=True).to_netcdf(path+savename)
        
        # Downwind SST gradient
        name='gradTheta_dwind'
        gradSST_dwind=xr.Dataset( {'gradTheta_dwind': ( ('time','lat'), np.zeros((lt,llat)) )},\
                coords={'time': dates,"lat": V0.lat.data } )
        for i in range(len(dates)):
                file=datetime.strftime(dates[i],folder2D+name+'_2D_'+'%Y%m%d_%H%M' + '.nc')
                gradSST_dwind.gradTheta_dwind[i,:]=xr.open_dataset(file).\
                        sel(lon=lonsel,method='nearest').sel(lat=slice(latmin,latmax))\
                                .gradTheta_dwind.data
        path=folderout+name+'_2Dv_'        
        gradSST_dwind.gradTheta_dwind.where(condi).mean(dim='time',skipna=True).to_netcdf(path+savename)
        
        # surface wind divergence      
        name='divU10'
        divU10=xr.Dataset( {'divU10': ( ('time','lat'), np.zeros((lt,llat)) )},\
                coords={'time': dates,"lat": V0.lat.data } )
        for i in range(len(dates)):
                file=datetime.strftime(dates[i],folder2D+name+'_2D_'+'%Y%m%d_%H%M' + '.nc')
                divU10.divU10[i,:]=xr.open_dataset(file).sel(lon=lonsel,method='nearest').sel(lat=slice(latmin,latmax)).divU10.data
                
        path=folderout+name+'_2Dv_'
        divU10.divU10.where(condi).mean(dim='time',skipna=True).to_netcdf(path+savename)

        # Vertical velocity averaged bellow 2 km
        name='W_2000mean'
        W_2000mean=xr.Dataset( {'W': ( ('time','lat'), np.zeros((lt,llat)) )},\
                coords={'time': dates,"lat": V0.lat.data } )
        for i in range(len(dates)):
                file=datetime.strftime(dates[i],folder2D+name+'_2D_'+'%Y%m%d_%H%M' + '.nc')
                W_2000mean.W[i,:]=xr.open_dataset(file).sel(lon=lonsel,method='nearest').sel(lat=slice(latmin,latmax)).W.data
                
        path=folderout+name+'_2Dv_'  
        W_2000mean.W.where(condi).mean(dim='time',skipna=True).to_netcdf(path+savename)
        
        
        
        ################# 3 - Extended Fig. 6
        
        # 3D Moisture flux convergence at a given longitude
        name='MFC'
        MFC=xr.Dataset( {'MFC': ( ('time','lat'), np.zeros((lt,llat)) )},\
                coords={'time': dates,"lat": V0.lat.data } )
        for i in range(len(dates)):
                file=datetime.strftime(dates[i],folder2D+name+'_2D_'+'%Y%m%d_%H%M' + '.nc')
                MFC.MFC[i,:]=xr.open_dataset(file).sel(lon=lonsel,method='nearest').sel(lat=slice(latmin,latmax)).MFC.data
        path=folderout+name+'_2Dv_'
        MFC.MFC.where(condi).mean(dim='time',skipna=True).to_netcdf(path+savename)
        
        # Total diabatic heating
        name='DTHDT_DELP'
        DTHDT_DELP=xr.Dataset( {'DTHDT_DELP': ( ('time','depth','lat'), np.zeros((lt,ldepth,llat)) )},\
                coords={'time': dates,"lat": V0.lat.data, "depth": np.array(list(V0.depth.data)) } )
        for i in range(len(dates)):
                file=datetime.strftime(dates[i],folder+name+'_2Dv_'+'%Y%m%d_%H%M' + '.nc')
                DTHDT_DELP.DTHDT_DELP[i,:]=xr.open_dataset(file).sel(lat=slice(latmin,latmax)).DTHDT_DELP.data
        path=folderout+name+'_2Dv_'
        DTHDT_DELP.DTHDT_DELP.where(condi).mean(dim='time',skipna=True).to_netcdf(path+savename)
        
        # Total vertical mass flux
        name='wrho'
        wrho=xr.Dataset( {'wrho': ( ('time','depth','lat'), np.zeros((lt,ldepth,llat)) )},\
                coords={'time': dates,"lat": V0.lat.data, "depth": np.array(list(V0.depth.data)) } )
        for i in range(len(dates)):
                file=datetime.strftime(dates[i],folder+name+'_2Dv_'+'%Y%m%d_%H%M' + '.nc')
                wrho.wrho[i,:]=xr.open_dataset(file).sel(lat=slice(latmin,latmax)).wrho.data
        path=folderout+name+'_2Dv_'
        wrho.wrho.where(condi).mean(dim='time',skipna=True).to_netcdf(path+savename)
        
        # Total precipitation
        name='PRECTOT'
        PRECTOT=xr.Dataset( {'PRECTOT': ( ('time','lat'), np.zeros((lt,llat)) )},\
                coords={'time': dates,"lat": V0.lat.data } )
        for i in range(len(dates)):
                file=datetime.strftime(dates[i],folder2D+name+'_2D_'+'%Y%m%d_%H%M' + '.nc')
                PRECTOT.PRECTOT[i,:]=xr.open_dataset(file).sel(lon=lonsel,method='nearest').sel(lat=slice(latmin,latmax)).PRECTOT.data
        path=folderout+name+'_2Dv_'  
        PRECTOT.PRECTOT.where(condi).mean(dim='time',skipna=True).to_netcdf(path+savename)


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
 
# Input data folder of 2D fields:      
folder2D='/path/to/your/extracted_data/2D/'

name,varname='V10M','V10M'
# getting shape
file=datetime.strftime(dates[0],folder2D+name+'_2D_'+'%Y%m%d_%H%M' + '.nc')
var0=xr.open_dataset(file)
lonn,latt=var0.lon.data,var0.lat.data
llat,llon=var0[varname].shape
lt=len(dates)
# creating dataset
var=xr.Dataset( {varname: ( ('time','lat','lon'), np.zeros((lt,llat,llon)) )},\
        coords={'time': dates,"lat": latt, "lon": lonn } )

# getting merdional wind values
for i in range(len(dates)):
        file=datetime.strftime(dates[i],folder2D+name+'_2D_'+'%Y%m%d_%H%M' + '.nc')
        var[varname][i,:]=xr.open_dataset(file)[varname].data
        
# mean meridional wind
V10M_mean=var[varname].mean(dim='lon').mean(dim='lat')


### 2 - Time mean with wind selection calculation

lonsel=152.5 # Meridional cross section

# Input folder of 3D fileds at a given longitude "2Dv" (meridional cross section)
# extracted with running_reading_GEOS_3D_lonsection.py
folder='/path/to/your/extracted_data/2Dv/'

# Output folder:
folderout='/path/to/your/extracted_data/2Dv/'

# Condition for wind selection
condi=(V10M_mean>=0) # Northward wind 
# condi=(V10M_mean<=0) # Southward wind

## Special handling for CNVMFC, which is stored at hh_30
shifted_dates=dates+timedelta(minutes=30)
V10M_mean_CNVMFC=V10M_mean.assign_coords(time= shifted_dates)
condi_CNVMFC=(V10M_mean_CNVMFC>=0) # Northward wind 
# condi_CNVMF=(V10M_mean_CNVMF<=0) # Southward wind

mean(dates,folder,folder2D,folderout,lonsel,condi,condi_CNVMFC,date_name+"_NWsup0_fullsection.nc")
# mean(dates,folder,lonsel,condi,"decjan_NWsup0_fullsection.nc")
# mean(dates,folder,lonsel,condi,"16dec21dec_SW_fullsection.nc")
# mean(dates,folder,lonsel,condi,"decjan_SW_fullsection.nc")

        
