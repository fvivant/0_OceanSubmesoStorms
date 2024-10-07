# Coded by Hector Torres (NASA-JPL)
# Adapted by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program for computing wavenumber spectra, cospectrum and coherence
# used in Extended Fig. 1
# Requirements: 
## 1) spectral_w.py

import numpy as np
import spactral_w as spactral_w
from scipy import signal
from numpy import pi
import time
from datetime import datetime,timedelta
import xarray as xr

#:::::::::::::::::::::::::::::::::::::::::::::::
def calc_ispec(k,l,E):
    """ calculates isotropic spectrum from 2D spectrum """

    dk,dl = k[1,]-k[0],l[1]-l[0]
    l,k = np.meshgrid(l,k)
    wv = np.sqrt(k**2 + l**2)

    if k.max()>l.max():
        kmax = l.max()
    else:
        kmax = k.max()

    # create radial wavenumber
    dkr = np.sqrt(dk**2 + dl**2)
    kr =  np.arange(dkr/2.,kmax+dkr,dkr)
    ispec = np.zeros(kr.size)

    for i in range(kr.size):
        fkr =  (wv>=kr[i]-dkr/2) & (wv<=kr[i]+dkr/2)
        # infinitesimal phase space angle in polar coordinates
        # dth = 2pi/fkr with fkr the discrete number of points 
        # in the disk bounded by two circle of radius kr-dkr/2 qnd kr+dkr/2
        # 2D case (x,y) of a w-omega spectrum with kx or ky not truncated for fft algorithm
        # fkr is the number over all the phase space
        # i.e. dth = pi / (fkr.sum())
        dth = pi / (fkr.sum())
        ispec[i] = E[fkr].sum() * kr[i] * dth *dkr

    return kr, ispec, dkr
#:::::::::::::::::::::::::::

# --------- Open files and turn them into numpy array ------------

dates = np.arange(datetime(2020,12,2,0),datetime(2021,2,28,0),timedelta(hours=1))\
       .astype(datetime)

# Input data path
path= '/path/to/your/extracted_data/2D/'

varnameA0="CAPE"
varnameA="CAPE"
varnameB0="Theta" #SST
varnameB="Theta"
filenameA=[varnameA0+'_2D_'+datetime.strftime(date,'%Y%m%d_%H%M')+'.nc' for date in dates]
filenameB=[varnameB0+'_3D_'+datetime.strftime(date,'%Y%m%d_%H%M')+'_ocean.nc' for date in dates]
fileref=xr.open_dataset(path+filenameA[0])
 
          
# ------- Symetric domain spectrum ---------------
nt,nlat,nlon=len(dates),fileref[varnameA].shape[0]*2,fileref[varnameA].shape[1]


A=np.zeros((nt,nlat,nlon))
B=np.zeros((nt,nlat,nlon))
t0=time.time()

for i in range(len(dates)):
    fileA=xr.open_dataset(path+filenameA[i])
    A[i,:nlat//2,:]=( np.array(fileA[varnameA].data)[:,:]) 
    A[i,nlat//2:,:]=( np.array(fileA[varnameA].data)[::-1,:])
    
    fileB=xr.open_dataset(path+filenameB[i])
    B[i,:nlat//2,:]=np.array(fileB[varnameB].data)[:,:]
    B[i,nlat//2:,:]=np.array(fileB[varnameB].data)[::-1,:]
print(time.time()-t0)


##-- This part is just to rearrange the array from time/lat/lon to lat/lon/time
A = np.swapaxes(A,0,2)
A = np.swapaxes(A,0,1)
#--

##-- This part is just to rearrange the array from time/lat/lon to lat/lon/time
B = np.swapaxes(B,0,2)
B = np.swapaxes(B,0,1)

# ------------ Detrending ----------------
A = signal.detrend(A,axis=0,type='linear')
A = signal.detrend(A,axis=1,type='linear')
B = signal.detrend(B,axis=0,type='linear')
B = signal.detrend(B,axis=1,type='linear')


iy,ix,it = A.shape
iaux =24
# iaux defines a window size to do statistics over spectra
# here iaux is a window of 1 day

nt = np.around(it/(iaux),decimals=1)
dx =np.around(0.07*np.pi*6371/180,decimals=2)

for i in range(int(nt)): # several iterations with these settings
    uauxA = A[:,:,i*iaux:i*iaux+iaux].mean(axis=2) # time mean of the window
    uauxB = B[:,:,i*iaux:i*iaux+iaux].mean(axis=2) # time mean of the window
    #i.e. uauxA and uauxB depend on space only (wavenumber spectra)

    if i == 0:
       EuA,EuB,cs,coh,k,l,d1,d2,Atofft,Btofft = spactral_w.spectra_w(uauxA,
                                                             uauxB,dx,dx)
    else:
       EuaA,EuaB,csa,coha,_,_,_,_,Atoffta,Btofftb = spactral_w.spectra_w(uauxA,
                                                             uauxB,dx,dx)
       EuA = EuA + EuaA
       EuB = EuB + EuaB
       cs = cs + csa
       coh = coh + coha
       Atofft=Atofft+Atoffta
       Btofft=Btofft+Btofftb

# The mean of spectra is considered here:
EuA = EuA/nt # Spectra of A
EuB = EuB/nt # Spectra of B
cs = cs/nt # Cospectrum of A and B
coh = coh/nt # Coherence of A and B
Atofft=Atofft/nt # A after detrending and hanning windowing
Btofft=Btofft/nt # B after detrending and hanning windowing

# ------------- Isospectra -----------------
kiso,EiA,dkiso = calc_ispec(k,l,EuA[:,:])
kiso,EiB,dkiso = calc_ispec(k,l,EuB[:,:])
kiso,Ei_cs,dkiso = calc_ispec(k,l,cs[:,:])
kiso,Ei_coh,dkiso = calc_ispec(k,l,coh[:,:])


# ---------- Save ouput --------------------
out= '/path/to/your/extracted_data/spectra_w/'
np.savez(out+'CAPE_SST_djf_d1day.npz',
         EiA=EiA,EiB=EiB,Ei_cs=Ei_cs,
         Ei_coh=Ei_coh,kiso=kiso)

print('------- End save --------')
print('------- Finish --------')

