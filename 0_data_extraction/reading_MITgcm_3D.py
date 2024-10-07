# Coded by Hector Torres (NASA-JPL)
# Adapted by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program required for extracting MITgcm data with running_3D_mitgcm.py
# Requirements: 
## 1) face_connections.py and llcmap_nea_split.py
## 2) "grid" folder available at:
## https://portal.nccs.nasa.gov/datashare/G5NR/DYAMONDv2/GEOS_6km_Atmosphere-MITgcm_4km_Ocean-Coupled/
## in the folder: MITgcm_output/grid
## 3) One must change "diro", "dirou", "dirov", line 188 and "GRIDDIR" by your own directories if you download the data
## 4) xgcm must be version 0.6.1 !!! Issues with 0.8.1

import numpy as np
import xarray as xr
from xmitgcm import open_mdsdataset
import xgcm
import warnings
warnings.filterwarnings("ignore")
from datetime import datetime,timedelta
from face_connections import face_connections
from llcmap_nea_split import LLCMap_nea_split


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


def mitgcm_xr(VAR,all_iters,lon,lat,level):

    dt = 45 # seconds 
        
    ### domain:
    dx = 0.07 ### GEOS interpolate onto 0.07 grid resolution.
    lat_out=np.arange(lat[0],lat[1]+dx,dx)
    lon_out=np.arange(lon[0],lon[1],dx)

    ##################################
    # MIT iterations
    ##################################

    iter0 = 0#120960
    delta_t=dt
    delta = 3600/delta_t # iters
    

    #### change this part
    diro = '/MITgcm_output/'+VAR+'/'
    ####

    ##################################
    # MIT grid info
    ##################################

    nx=2160
    #### change this part
    GRIDDIR='/MITgcm_output/grid/'
    ####

    ds = open_mdsdataset(diro, grid_dir=GRIDDIR, iters=all_iters, 
                       geometry='llc', read_grid=True, 
                       default_dtype=np.dtype('>f4'),
                       delta_t=delta_t, ignore_unknown_vars=True, nx=nx)

    grid = xgcm.Grid(ds, periodic=False, face_connections=face_connections)

    coords = ds.coords.to_dataset().reset_coords()
    msk=coords.hFacC.sel(k=0)
    msk=msk.where(msk>0,np.nan) 
    XC=coords.XC*msk
    YC=coords.YC*msk
    mapper = LLCMap_nea_split(YC.values, XC.values,lat_out,
                            lon_out,radius=10e3)


    ds = open_mdsdataset(diro, 
        grid_dir=GRIDDIR, 
                         iters=all_iters, geometry='llc',
                         read_grid=True, 
                         default_dtype=np.dtype('>f4'), 
                         delta_t=delta_t, 
                         ignore_unknown_vars=True, nx=nx)
    

    if VAR == 'V':
        #### change this part
        dirou='/MITgcm_output/U/' #U is required to extract V

        ####
        dsu = open_mdsdataset(dirou,
                      grid_dir=GRIDDIR,
                            iters=all_iters,geometry='llc',read_grid=True,
                            default_dtype=np.dtype('>f4'),delta_t=delta_t,
                            ignore_unknown_vars=True,nx=nx)
        
        ds1 = open_mdsdataset(GRIDDIR, iters=1,
                            ignore_unknown_vars=True, geometry='llc', nx=nx)

        AngleCS=ds1['CS'];AngleSN=ds1['SN'];
        UV=grid.interp_2d_vector({'X': dsu['U'].isel(k=level), 
                    'Y': ds['V'].isel(k=level)},boundary='fill')
        x=(UV['X']*AngleSN+UV['Y']*AngleCS)


    elif VAR == 'U':
        #### change this part
        dirov='/MITgcm_output/V/'
        ####
        dsv = open_mdsdataset(dirov,
                      grid_dir=GRIDDIR,
                            iters=all_iters,geometry='llc',read_grid=True,
                            default_dtype=np.dtype('>f4'),delta_t=delta_t,
                            ignore_unknown_vars=True,nx=nx)

        ds1 = open_mdsdataset(GRIDDIR, iters=1, 
                            ignore_unknown_vars=True, 
                            geometry='llc', nx=nx) 
        AngleCS=ds1['CS'];AngleSN=ds1['SN'];
        UV=grid.interp_2d_vector({'X': ds['U'].isel(k=level), 'Y': dsv['V'].isel(k=level)},boundary='fill')
        x=(UV['X']*AngleCS-UV['Y']*AngleSN)

    elif VAR == 'Theta':
        x = ds['Theta'].isel(k=level)

    elif VAR == 'W':
        x=grid.interp(ds['W'],'Z',to='center',boundary='fill').isel(k=level)

    elif VAR == 'Eta':
        x = ds[VAR]

    ##################################################
    ##### Mapping X to new grid 
    TMP=xr.DataArray(mapper(x.values))
    
    return TMP,lon_out,lat_out

def getting_3D(VAR,all_iters,LON,LAT,levels):
    import scipy.io
    ## reading z-coordinates
    #### change this part
    GRIDDIR='/MITgcm_output/grid/'
    ####
    z = scipy.io.loadmat(GRIDDIR+'thk90.mat')
    dpt = z['dpt90'][...]
    zaxis = np.empty((len(levels)))
    I = 0

    for i in range(0,len(levels)):
        zaxis[i] = dpt[:,i]
        print(I)
        data,lon,lat = mitgcm_xr(VAR,all_iters,LON,LAT,levels[i])
        if I == 0:
            MAT = np.empty((len(levels),len(lat),len(lon)))

        MAT[I,:,:] = data
        I += 1
    return MAT,lon,lat,zaxis
