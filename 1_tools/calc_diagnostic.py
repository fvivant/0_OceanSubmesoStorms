# Coded by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program for computing diagnostic variables such as divU10, gradSST, moisture flux convergence

import xarray as xr
import numpy as np
from datetime import datetime,timedelta

#-----------------------------------------------
y1,m1,d1,h1,M1=2020,12,1,0,0
y2,m2,d2,h2,M2=2021,2,28,0,0

pathin='/path/to/your/extracted_data/2D/'
pathout='/path/to/your/extracted_data/2D/'

dH=1
date = np.arange(datetime(y1,m1,d1,h1,M1),datetime(y2,m2,d2,h2,M2),timedelta(hours=dH))\
       .astype(datetime)

Rt=6356.7523e3 #earthradius

for i in range(len(date)):
       ## SST gradient:
       sst=xr.open_dataset(pathout+"Theta_3D_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'_ocean.nc')
       gradx,grady=sst.Theta.differentiate("lon").data*180/(np.pi*Rt),sst.Theta.differentiate("lat").data*180/(np.pi*Rt)
       gradsst=xr.Dataset( { "gradTheta": ( ('lat','lon'), (gradx**2+grady**2)**0.5 ) }, coords={'lat': sst.lat.data, 'lon': sst.lon.data}, )
       gradsst.to_netcdf(pathout+"gradTheta_2D_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'_ocean.nc')

       # ## Wind divergence at 10m:
       # u10=xr.open_dataset(pathout+"U10M_2D_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc')
       # v10=xr.open_dataset(pathout+"V10M_2D_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc')
       # dxU,dyV=u10.U10M.differentiate("lon").data*180/(np.pi*Rt),v10.V10M.differentiate("lat").data*180/(np.pi*Rt)
       # divU10=xr.Dataset( { "divU10": ( ('lat','lon'), dxU+dyV ) }, coords={'lat': u10.lat.data, 'lon': u10.lon.data}, )
       # divU10.to_netcdf(pathout+"divU10_2D_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc')

       # ## Downwind SST gradient:
       # sst=xr.open_dataset(pathout+"Theta_3D_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'_ocean.nc')
       # gradx,grady=sst.Theta.differentiate("lon").data*180/(np.pi*Rt),sst.Theta.differentiate("lat").data*180/(np.pi*Rt)
       # u10=xr.open_dataset(pathout+"U10M_2D_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc')
       # v10=xr.open_dataset(pathout+"V10M_2D_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc')
       # gradTheta_dwind=xr.Dataset( { "gradTheta_dwind": ( ('lat','lon'),\
       #               (gradx*u10.U10M.data+grady*v10.V10M.data)/(u10.U10M.data**2+v10.V10M.data**2)**0.5 ) },\
       #                      coords={'lat': u10.lat.data, 'lon': u10.lon.data}, )
       # gradTheta_dwind.to_netcdf(pathout+"gradTheta_dwind_2D_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc')

       # ## Moisture flux convergence at 10m (specific humidity at 2 meters is used):
       # QV2M=xr.open_dataset(pathout+"QV2M_2D_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc')
       # divU10=xr.open_dataset(pathout+"divU10_2D_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc')
       # u10=xr.open_dataset(pathout+"U10M_2D_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc')
       # v10=xr.open_dataset(pathout+"V10M_2D_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc')
       # dxq,dyq=QV2M.QV2M.differentiate("lon").data*180/(np.pi*Rt),QV2M.QV2M.differentiate("lat").data*180/(np.pi*Rt)
       # mfc=-u10.U10M*dxq-v10.V10M*dyq -QV2M.QV2M*divU10.divU10
       # MFC= xr.Dataset( { "MFC": ( ('lat','lon'),  mfc.data) }, coords={'lat': QV2M.lat.data, 'lon': QV2M.lon.data}, )
       # MFC.to_netcdf(pathout+"MFC_2D_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc')

       # ## MITgcm vorticity
       # U=xr.open_dataset(pathout+"U_3D_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'_ocean.nc')
       # V=xr.open_dataset(pathout+"V_3D_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'_ocean.nc')
       # udy,vdx=U.U.differentiate("lat").data*180/(np.pi*Rt),V.V.differentiate("lon").data*180/(np.pi*Rt)
       # vort=xr.Dataset( { "vort": ( ('lat','lon'), vdx-udy ) }, coords={'lat': U.lat.data, 'lon': U.lon.data}, )
       # vort.to_netcdf(pathout+"vort_2D_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'_ocean.nc')

       # ### Density and total convectif mass flux:
       # T=xr.open_dataset(pathout+"T_2Dv_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc')
       # P=xr.open_dataset(pathout+"P_2Dv_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc')
       # rho= xr.Dataset( { "rho": ( ('depth','lat'),  P.P.data/(287*T.T.data)) }, coords={'depth': T.depth.data,'lat': T.lat.data})
       # rho.to_netcdf(pathout+"rho_2Dv_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc')
       # W=xr.open_dataset(pathout+"W_2Dv_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc')
       # wrho= xr.Dataset( { "wrho": ( ('depth','lat'),  rho.rho.data*W.W.data) }, coords={'depth': T.depth.data,'lat': T.lat.data})
       # wrho.to_netcdf(pathout+"wrho_2Dv_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc')

       # ## DTHDT (total diabatic heating) de-weighted with pressure thickness:
       # DTHDT=xr.open_dataset(pathout+"DTHDT_2Dv_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc')
       # DELP=xr.open_dataset(pathout+"DELP_2Dv_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc')
       # DTHDT_DELP=xr.Dataset( { "DTHDT_DELP": ( ('depth','lat'), DTHDT.DTHDT.data/DELP.DELP.data ) }, coords={ 'depth': DTHDT.depth.data,'lat': DTHDT.lat.data}, )
       # DTHDT_DELP.to_netcdf(pathout+"DTHDT_DELP_2Dv_"+datetime.strftime(date[i],'%Y%m%d_%H%M')+'.nc')
       


