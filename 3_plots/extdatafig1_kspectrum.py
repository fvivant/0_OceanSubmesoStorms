# Coded by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program to plot Extended data Fig.1

import matplotlib.pyplot as plt
import numpy as np
from numpy import load


out = '/path/to/your/extracted_data/spectra_w/'
data = load(out+'CAPE_SST_djf_d1h.npz')
lst = data.files
EiA=data['EiA']
EiB=data['EiB']
Ei_cs=data['Ei_cs']
Ei_coh=data['Ei_coh']
kiso=data['kiso']


#PLOT
size=30
parameters = {'axes.labelsize': size, 'axes.titlesize': size,"axes.labelsize" : size,"xtick.labelsize" : size,\
       "ytick.labelsize" : size}
plt.rcParams.update(parameters)


fig=plt.figure(figsize=(22,15))
fig.subplots_adjust(bottom=0.09, top=0.87, left=0.09, right=0.92,
                    wspace=0.05, hspace=0.25)
axx= plt.subplot(111)
axx2 = axx.twinx()
        
lambdadisp=np.array([20,50,70,100,200,500,700,1000])
kdisp=1/lambdadisp
ymin,ymax=0,1.15#1e0
xmin,xmax=4e-4,kiso.max()
lin=axx2.vlines(x=kdisp, ymin=ymin, ymax=ymax, color = 'k',linestyle='dashed',\
    linewidth=2)
lin.set_clip_on(False)
for kk in range(len(kdisp)):
    lin=axx2.text(kdisp[kk], 1.01, str(lambdadisp[kk])+'km',color='k',\
        rotation=90, verticalalignment='bottom',horizontalalignment='right',size=size)
    lin.set_clip_on(False)
axx2.set_ylim(0,1)
# --------------------------
lw=5

axx.loglog(kiso,(EiA*kiso),c='r',linewidth=lw,label='CAPE snapshot')
Ei_coh2=Ei_cs/(np.sqrt(EiA)*np.sqrt(EiB))
axx2.semilogx(kiso,(Ei_coh2),c='r',linewidth=lw,linestyle='dashed',label='Coherence snapshot')#/(Ei_cs*kiso).max()

# --------------------------
out = '/path/to/your/extracted_data/spectra_wf/'
data = load(out+'CAPE_SST_djf_d1day.npz')
lst = data.files
EiA=data['EiA']
EiB=data['EiB']
Ei_cs=data['Ei_cs']
Ei_coh=data['Ei_coh']
kiso=data['kiso']


axx.loglog(kiso,(EiA*kiso),c='b',linewidth=lw,label='CAPE 1-day mean')#linestyle='dashed'

Ei_coh2=Ei_cs/(np.sqrt(EiA)*np.sqrt(EiB))
axx2.semilogx(kiso,(Ei_coh2),c='b',linestyle='dashed',linewidth=lw,label='Coherence 1-day mean')

#--------------------
axx.set_ylabel('$\mathcal{I}(|\widehat{CAPE}|^{2}) \\times |\mathbf{k}|$ [(J.kg$^{-1}$.km)$^{2}$.cpkm]',labelpad=10)
axx2.set_ylabel('Coherence with SST',labelpad=20)    
axx.set_xlabel('$|\mathbf{k}|$ [cpkm]') 

axx2.legend(loc='upper right', prop={'size': 25})
axx.legend(loc='upper left', prop={'size': 25})

plt.savefig('/path/to/your/plots/ExtDataFig1.png')

plt.close()


