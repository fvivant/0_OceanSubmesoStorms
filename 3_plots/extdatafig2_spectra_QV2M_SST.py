# Coded by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program to plot Extended data Fig.2

import numpy as np
import pylab as plt
import matplotlib.colors as mcolors
from numpy import load

out = '/path/to/your/extracted_data/spectra_wf/'
data = load(out+'QV2M_SST_symdom.npz')
lst = data.files
EisoA=data['EisoA']
EisoB=data['EisoB']
EisoCS=data['EisoCS']
EisoCOH=data['EisoCOH']
kiso=data['kiso']
om=data['om']
Atofft=data['Atofft']
Btofft=data['Btofft']
dx=data['dx']
dt=data['dt']

def Parseval_coef(spec,om,kiso):
    perseval_isospec=np.trapz(np.sum(spec,axis=0),x=om,axis=0)*2
    
    Om,Kiso=np.meshgrid(om,kiso)
    period=1/(Om*24)
    # Phase space segmentation
    #Synoptic
    condi=np.logical_and(1/Kiso>=200,period<5)
    pers1=np.trapz(np.sum(np.where(condi,spec,spec*0),axis=0),x=om,axis=0)*2

    #Meso
    condi1=np.logical_and(1/Kiso>=70,1/Kiso<=700)
    condi=np.logical_and(condi1,period>=5)
    pers2=np.trapz(np.sum(np.where(condi,spec,spec*0),axis=0),x=om,axis=0)*2

    #Submeso
    condi=np.logical_and(1/Kiso<70,period>=1)
    pers3=np.trapz(np.sum(np.where(condi,spec,spec*0),axis=0),x=om,axis=0)*2

    #HighAtm
    condi=np.logical_and(1/Kiso<200,period<1)
    pers4=np.trapz(np.sum(np.where(condi,spec,spec*0),axis=0),x=om,axis=0)*2

    #Overlap
    condi1=np.logical_and(1/Kiso>70,1/Kiso<200)
    condi2=np.logical_and(period>1,period<5)
    condi=np.logical_and(condi1,condi2)
    pers5=np.trapz(np.sum(np.where(condi,spec,spec*0),axis=0),x=om,axis=0)*2
    
    #Ocean large
    condi=np.logical_and(1/Kiso>700,period>5)
    pers6=np.trapz(np.sum(np.where(condi,spec,spec*0),axis=0),x=om,axis=0)*2
    return perseval_isospec,pers1,pers2,pers3,pers4,pers5,pers6

def plot_Perseval_coef(ax,perseval_isospec,pers1,pers2,pers3,pers4,pers5,pers6,kdisp,ymin,ymax,xmin,xmax):
    size=6
    pad=0.1

    ### horizontal
    # Synoptic
    ax.hlines(y=1/(5*24), xmin=xmin, xmax=1/70, color = 'k',linestyle='dashed',\
        linewidth=size)
    ax.text( 1/200 *(1-pad),1/(5*24) *(1+pad), str(round(100*pers1/perseval_isospec))+'%',color='k',
            rotation=0, verticalalignment='bottom',horizontalalignment='right',size=40)
    #Highatm
    ax.hlines(y=1/(1*24), xmin=1/200, xmax=xmax, color = 'k',linestyle='dashed',\
        linewidth=size)
    ax.text( 1/200 *(1+pad),1/(1*24) *(1+pad), str(round(pers4*100/perseval_isospec))+'%',color='k',rotation=0, verticalalignment='bottom',horizontalalignment='left',size=40)

    ### Vertical
    #Synoptic
    ax.vlines(x=1/200, ymin=1/(5*24), ymax=ymax, color = 'k',linestyle='dashed',\
        linewidth=size)
    # ax.text( 1/100 *(1+pad),1/(5*24)*(1-pad), str(int(pers3*100/perseval_isospec))+'%',color='k',rotation=0, verticalalignment='top',horizontalalignment='left',size=40)

    ax.vlines(x=1/70, ymin=ymin, ymax=1/(1*24), color = 'k',linestyle='dashed',\
        linewidth=size)
    #meso
    ax.text( 1/70 *(1-pad),1/(5*24)*(1-pad), str(round(pers2*100/perseval_isospec))+'%',color='k',rotation=0, verticalalignment='top',horizontalalignment='right',size=40)
    #submeso
    ax.text( 1/70 *(1+pad),1/(1*24)*(1-pad), str(round(pers3*100/perseval_isospec))+'%',color='k',rotation=0, verticalalignment='top',horizontalalignment='left',size=40)

    #Ocean large
    ax.vlines(x=1/700, ymin=ymin, ymax=1/(5*24), color = 'k',linestyle='dashed',\
        linewidth=size)
    ax.text( 1/1000 *(1-pad),1/(20*24)*(1-pad), str(round(pers6*100/perseval_isospec))+'%',color='k',rotation=0, verticalalignment='top',horizontalalignment='right',size=35)
    
    #Overlap
    ax.text( 1/100,1/(2.5*24), str(round(pers5*100/perseval_isospec))+'%',color='k',rotation=0, verticalalignment='bottom',horizontalalignment='right',size=35)


    
size=30
parameters = {'axes.labelsize': size, 'axes.titlesize': size,"axes.labelsize" : size,"xtick.labelsize" : size,\
       "ytick.labelsize" : size,"ytick.major.size" : 8,"ytick.major.width" : 3,"ytick.minor.size" : 4,"ytick.minor.width" : 2,\
       "xtick.major.size" : 8,"xtick.major.width" : 3,"xtick.minor.size" : 4,"xtick.minor.width" : 2}
plt.rcParams.update(parameters)
nrows,ncols=1,3
fig, axes=plt.subplots(nrows,ncols,sharey=True,sharex=True,figsize=(40,15))#,constrained_layout=True)

    
fig.subplots_adjust(bottom=0.1, top=0.85, left=0.06, right=0.98,
                    wspace=0.05, hspace=0)

tx,ty=0.04,0.95
#------------------------------------------------------------------------------ 
i,j=0,0
axes[j].text(tx,ty,'a',weight="bold",size=35,horizontalalignment='center',verticalalignment='center',\
                transform = axes[j].transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
omlim=1
ymin,ymax=4e-4,om.max()#1e0
xmin,xmax=4e-4,kiso.max()
lambdadisp=np.array([20,50,100,200,500,1000])
kdisp=1/lambdadisp
Tdisp=np.array([1,5,10,30])
omdisp=1/(Tdisp*24)

conserv=kiso[None,:]*om[om<omlim][:,None]

disp=EisoA.T[:,:]*conserv
specmin=disp.max()*1e-3
# print(disp.max())
norm=mcolors.LogNorm(vmin=specmin)
set=axes[j].pcolormesh(kiso,om,disp,\
    cmap='jet',norm=norm,alpha=0.85)

axes[j].vlines(x=kdisp, ymin=ymin, ymax=ymax, color = 'k',linestyle='dashed',\
    linewidth=2)
for kk in range(len(kdisp)):
    axes[j].text(kdisp[kk], ymax, str(lambdadisp[kk])+'km',color='k',\
        rotation=90, verticalalignment='top',horizontalalignment='right',size=30)
axes[j].hlines(y=omdisp, xmin=xmin, xmax=xmax, color = 'k',linestyle='dashed',\
    linewidth=2)
for kk in range(len(Tdisp)):
    axes[j].text( xmax,omdisp[kk], str(Tdisp[kk])+'days',color='k',\
        rotation=0, verticalalignment='bottom',horizontalalignment='right',size=30)
    
dcolor=0.015
tick_font_size=30
labelpad=20
pos1 = axes[j].get_position()
cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
cbar=plt.colorbar(mappable=set,cax=cax,orientation="horizontal",ticklocation='top',aspect=35,extend='max')
cbar.set_label("$\mathcal{I}(|\widehat{q_{2m}}|^{2}) \\times f \\times |\mathbf{k}|$ [(kg.kg$^{-1}$.km.h)$^{2}$.cph.cpkm]",labelpad=labelpad)
cbar.ax.tick_params(labelsize=tick_font_size)

axes[j].set_ylim(ymin,ymax)
axes[j].set_xlim(xmin,xmax)
axes[j].set_ylabel('$f$ [cph]')
axes[j].set_xlabel('$|\mathbf{k}|$ [cpkm]')
# axes[i,j].set_xlabel('k')
axes[j].set_xscale('log')
axes[j].set_yscale('log')

perseval_isospec,pers1,pers2,pers3,pers4,pers5,pers6=Parseval_coef(EisoA,om,kiso)
plot_Perseval_coef(axes[j],perseval_isospec,pers1,pers2,pers3,pers4,pers5,pers6,kdisp,ymin,ymax,xmin,xmax)

#------------------------------------------------------------------------------ 
i,j=0,1

axes[j].text(tx,ty,'b',weight="bold",size=35,horizontalalignment='center',verticalalignment='center',\
                transform = axes[j].transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))

disp=EisoB.T[:,:]*conserv
specmin=disp.max()*1e-3
print(disp.max())
norm=mcolors.LogNorm(vmin=specmin)
set=axes[j].pcolormesh(kiso,om[om<omlim],disp
                       ,cmap='jet',norm=norm,alpha=0.85)#np.linspace(0,5e-3,30))
axes[j].vlines(x=kdisp, ymin=ymin, ymax=ymax, color = 'k',linestyle='dashed',\
    linewidth=2)
for kk in range(len(kdisp)):
    axes[j].text(kdisp[kk], ymax, str(lambdadisp[kk])+'km',color='k',\
        rotation=90, verticalalignment='top',horizontalalignment='right',size=30)
axes[j].hlines(y=omdisp, xmin=xmin, xmax=xmax, color = 'k',linestyle='dashed',\
    linewidth=2)
for kk in range(len(Tdisp)):
    axes[j].text( xmax,omdisp[kk], str(Tdisp[kk])+'days',color='k',\
        rotation=0, verticalalignment='bottom',horizontalalignment='right',size=30)

pos1 = axes[j].get_position()
cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
cbar=plt.colorbar(mappable=set,cax=cax,orientation="horizontal",ticklocation='top',aspect=35,extend='max')
cbar.set_label("$\mathcal{I}(|\widehat{SST}|^{2}) \\times f \\times |\mathbf{k}|$ [(°C.km.h)$^{2}$.cph.cpkm]",labelpad=labelpad)
cbar.ax.tick_params(labelsize=tick_font_size)

axes[j].set_xlabel('$|\mathbf{k}|$ [cpkm]')

perseval_isospec,pers1,pers2,pers3,pers4,pers5,pers6=Parseval_coef(EisoB,om,kiso)
print((pers1+pers2+pers3+pers4+pers5+pers6)/perseval_isospec)
plot_Perseval_coef(axes[j],perseval_isospec,pers1,pers2,pers3,pers4,pers5,pers6,kdisp,ymin,ymax,xmin,xmax)


#------------------------------------------------------------------------------ 
i,j=0,2

axes[j].text(tx,ty,'c',weight="bold",size=35,horizontalalignment='center',verticalalignment='center',\
                transform = axes[j].transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))

field=EisoCS.T[om<omlim,:]*conserv
norm=mcolors.SymLogNorm(linthresh=field.max()*1e-1,
                                vmin=-field.max(),vmax=field.max(), base=10)

set=axes[j].pcolormesh(kiso,om[om<omlim],field
                       ,cmap='seismic',norm=norm)#np.linspace(0,5e-3,30))
axes[j].vlines(x=kdisp, ymin=ymin, ymax=ymax, color = 'k',linestyle='dashed',\
    linewidth=2)
for kk in range(len(kdisp)):
    axes[j].text(kdisp[kk], ymax, str(lambdadisp[kk])+'km',color='k',\
        rotation=90, verticalalignment='top',horizontalalignment='right',size=30)
axes[j].hlines(y=omdisp, xmin=xmin, xmax=xmax, color = 'k',linestyle='dashed',\
    linewidth=2)
for kk in range(len(Tdisp)):
    axes[j].text( xmax,omdisp[kk], str(Tdisp[kk])+'days',color='k',\
        rotation=0, verticalalignment='bottom',horizontalalignment='right',size=30)

pos1 = axes[j].get_position()
cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
cbar=plt.colorbar(mappable=set,cax=cax,orientation="horizontal",ticklocation='top',aspect=35,extend='both')
cbar.set_label("$\mathcal{I}(\widehat{q_{2m}}.\widehat{SST}^*) \\times f \\times |\mathbf{k}|$ [kg.kg$^{-1}$.°C.(km.h)$^{2}$.cph.cpkm]",labelpad=labelpad)
cbar.ax.tick_params(labelsize=tick_font_size)

axes[j].set_xlabel('$|\mathbf{k}|$ [cpkm]')

perseval_isospec,pers1,pers2,pers3,pers4,pers5,pers6=Parseval_coef(EisoCS,om,kiso)
# print((pers1+pers2+pers3+pers4+pers5+pers6)/perseval_isospec)
plot_Perseval_coef(axes[j],perseval_isospec,pers1,pers2,pers3,pers4,pers5,pers6,kdisp,ymin,ymax,xmin,xmax)

plt.savefig('/path/to/your/plots/ExtDataFig2.png')

