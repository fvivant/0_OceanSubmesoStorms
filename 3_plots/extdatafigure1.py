# Coded by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program to plot Extended data Fig.1

import matplotlib.pyplot as plt
import numpy as np
from numpy import load


# out = '/path/to/your/extracted_data/spectra_w/'
out = '/nobackup/fvivant/DATA/KExt/winter_decjanfev/spectrum/'
data = load(out+'CAPE_SST_djf_d1h.npz')
lst = data.files
EiA=data['EiA']
EiB=data['EiB']
Ei_cs=data['Ei_cs']
Ei_coh=data['Ei_coh']
kiso=data['kiso']


#PLOT
size=30
parameters = {'axes.labelsize': size, 'axes.titlesize': size,"axes.labelsize" : size,"xtick.labelsize" : size,
       "ytick.labelsize" : size}
plt.rcParams.update(parameters)

figsize = (22,15)
fig=plt.figure(figsize=figsize)
fig.subplots_adjust(bottom=0.09, top=0.87, left=0.09, right=0.92,
                    wspace=0.05, hspace=0.25)
axx= plt.subplot(111)
axx2 = axx.twinx()
        
lambdadisp=np.array([20,50,70,100,200,500,700,1000])
kdisp=1/lambdadisp
ymin,ymax=0,1.15#1e0
xmin,xmax=4e-4,kiso.max()
lin=axx2.vlines(x=kdisp, ymin=ymin, ymax=ymax, color = 'k',linestyle='dashed',
    linewidth=2)
lin.set_clip_on(False)
for kk in range(len(kdisp)):
    lin=axx2.text(kdisp[kk], 1.01, str(lambdadisp[kk])+'km',color='k',
        rotation=90, verticalalignment='bottom',horizontalalignment='right',size=size)
    lin.set_clip_on(False)
axx2.set_ylim(0,1)
# --------------------------
lw=5

axx.loglog(kiso,(EiA*kiso),c='r',linewidth=lw,label='CAPE snapshot')
Ei_coh2=Ei_cs/(np.sqrt(EiA)*np.sqrt(EiB))
axx2.semilogx(kiso,(Ei_coh2),c='r',linewidth=lw,linestyle='dashed',label='Coherence snapshot')#/(Ei_cs*kiso).max()

# --------------------------
# out = '/path/to/your/extracted_data/spectra_wf/'
out = '/nobackup/fvivant/DATA/KExt/winter_decjanfev/spectrum/'
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
axx.set_ylabel('$\mathcal{I}(|\widehat{\mathrm{CAPE}}|^{2}) \\times |\mathbf{k}|$ [(J.kg$^{-1}$.km)$^{2}$.cpkm]',labelpad=10)
axx2.set_ylabel('Coherence with SST',labelpad=20)    
axx.set_xlabel('$|\mathbf{k}|$ [cpkm]') 

axx2.legend(loc='upper right', prop={'size': 25})
axx.legend(loc='upper left', prop={'size': 25})


fig_width = figsize[0]  # inches in this program
fig_width_required = 18  # cm required by the journal
cm = 2.54  # inches to cm
scaling = cm * fig_width / fig_width_required

dpi_target = 500  # desired dpi
dpi = dpi_target / scaling  # to get a target dpi after resizing
figname = 'ExtDataFigure1.png'
path = '/nobackup/fvivant/programs/1_github/Figures/'
# path='/path/to/your/plots/'

plt.savefig('/nobackup/fvivant/programs/1_github/Figures/' + figname,
            format=figname[-3:],
            dpi=dpi)
plt.close()

resize_plot = True
# ----------------- Rescaling ------------------------------------------
if resize_plot == True and figname[-3:] == 'png':

        from PIL import Image

        # Open the image
        input_file = path + figname
        output_file = path + figname

        # Open the image
        image = Image.open(input_file)

        # Save with new DPI
        image.save(output_file, dpi=(dpi_target, dpi_target))

        # Get pixel dimensions
        image = Image.open(output_file)
        width, height = image.size
        # Get DPI (if available)
        dpi = image.info.get("dpi")
        print(
            f"Width: {width*cm/dpi[0]} cm, Height: {height*cm/dpi[1]} cm, DPI: {dpi[0]} x {dpi[1]}"
        )


