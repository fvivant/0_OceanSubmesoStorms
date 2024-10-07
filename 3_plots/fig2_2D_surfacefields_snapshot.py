# Coded by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program to plot Fig.2

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import matplotlib.gridspec as gridspec
from metpy.plots import ColdFront
import cmocean

# ------------------------------ Plot function ---------------------
def plott(endname):
        # Function to plot Figure 2
        # Input: endname is the common name at the end of each variable file 
        # e.g. '_2D_20201227_1300.nc'
        
        # ------------- Open data
        lonmin,lonmax=130,170
        latmin,latmax=25,50
        
        # larger domain than the domain of study for panel a
        path='/path/to/your/extracted_data/2D_large/'
        
        name='U10M'
        file=path+name+endname
        U10Ml=xr.open_dataset(file).sel(lon=slice(lonmin,lonmax)).sel(lat=slice(latmin,latmax))
        
        name='V10M'
        file=path+name+endname
        V10Ml=xr.open_dataset(file).sel(lon=slice(lonmin,lonmax)).sel(lat=slice(latmin,latmax))
        
        #SST from MITgcm
        name='Theta'
        file=path+name+'_3'+endname[2:-3]+'_ocean.nc'
        SSTl=xr.open_dataset(file).isel(depth=0).sel(lon=slice(lonmin,lonmax)).sel(lat=slice(latmin,latmax))

        # domain of study for panel b, c, d
        path='/path/to/your/extracted_data/2D/'
        
        name='Theta'
        file=path+name+'_3'+endname[2:-3]+'_ocean.nc'
        SST2=xr.open_dataset(file)
        
        name='vort'
        file=path+name+endname[:-3]+'_ocean.nc'
        vort=xr.open_dataset(file)

        name='U10M'
        file=path+name+endname
        U10M=xr.open_dataset(file)
        
        name='V10M'
        file=path+name+endname
        V10M=xr.open_dataset(file)
        
        name='EFLUX'
        file=path+name+endname
        EFLUX=xr.open_dataset(file)
        
        name='divU10'
        file=path+name+endname
        divU10=xr.open_dataset(file)

        # -------------- Definition of the plot        
        size=20
        parameters = {'axes.labelsize': size, 'axes.titlesize': size,"axes.labelsize" : size,"xtick.labelsize" : size,\
                "ytick.labelsize" : size,"ytick.major.size" : 5,"ytick.major.width" : 2,"ytick.minor.size" : 3,"ytick.minor.width" : 1,\
       "xtick.major.size" : 5,"xtick.major.width" : 2,"xtick.minor.size" : 3,"xtick.minor.width" : 1}
        plt.rcParams.update(parameters)
        
        fig=plt.figure(figsize=(20,20))
        
        gs=gridspec.GridSpec(2,2)
        gs.update(bottom=0.03, top=0.93, left=0.05, right=0.98,
                        wspace=0.12, hspace=0.25)

        # --------------------------------------------------------------
        axx=plt.subplot(gs[0,0])
        
        nlev=30
        tx,ty=0.04,0.95
        axx.text(tx,ty,'a',weight="bold",size=25,horizontalalignment='center',verticalalignment='center',\
                transform = axx.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        axx.text(138.8,36.4,'Japan',weight="bold",size=20,horizontalalignment='center',verticalalignment='center',rotation=40)

        #Cold front
        angle=50
        axx.plot([154.2,148.7],[31.2+np.cos(angle)*(154.2-148.7),31.2],c='b',linewidth=3,path_effects=[ColdFront(size=6, spacing=1.5)])
        #SST
        c2=(SSTl.Theta).plot.contourf(ax=axx,\
                levels=np.arange(0,26,0.5),cmap='coolwarm',alpha=1,add_colorbar=False,extend='max')
        
        #Wind field
        lonnn,lattt=np.meshgrid(U10Ml.lon.data,U10Ml.lat.data)             
        lkey=20
        alpha=15
        quiverscale=None
        quiverwidth=0.003
        
        q=axx.quiver(lonnn[::alpha,::alpha],lattt[::alpha,::alpha],\
                U10Ml.U10M.data[::alpha,::alpha],V10Ml.V10M.data[::alpha,::alpha],zorder=2,\
                width=quiverwidth,scale=quiverscale,angles='xy', scale_units='xy')
        ##wind field legend
        txq,tyq=0.86,0.96
        tt=axx.text(txq+0.04,tyq,'               ',weight="bold",size=20,horizontalalignment='center',verticalalignment='center',\
                transform = axx.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        tt.set_zorder(2)
        qt=axx.quiverkey(q, X=txq, Y=tyq, U=lkey,
             label=str(lkey)+' m.s$^{-1}$', labelpos='E',fontproperties={'weight': 'bold','size': 15},zorder=3)#,bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        t = qt.text.set_zorder(3)
        
        #Domain of study
        loni=[142.3,162]
        lati=[31,41.85]
        co='blue'
        linestyle='solid'
        linewidth=4
        axx.vlines(loni[0],ymin=lati[0],ymax=lati[1],colors=co,linestyles=linestyle,linewidths=linewidth)
        axx.vlines(loni[1],ymin=lati[0],ymax=lati[1],colors=co,linestyles=linestyle,linewidths=linewidth)
        axx.hlines(lati[0],xmin=loni[0],xmax=loni[1],colors=co,linestyles=linestyle,linewidths=linewidth)
        axx.hlines(lati[1],xmin=loni[0],xmax=loni[1],colors=co,linestyles=linestyle,linewidths=linewidth)
        
        #Typical size
        axx.annotate('',size=15,xy=(loni[1], lati[1]), xytext=(175, 50),ha='center', va='top',
            arrowprops=dict(arrowstyle='-',color=co, lw=linewidth))
        axx.annotate('',size=15,xy=(loni[1], lati[0]), xytext=(175, 25),ha='center', va='top',
            arrowprops=dict(arrowstyle='-',color=co, lw=linewidth))
        
        # Colorbar, ticks and labels
        tick_font_size=18
        labelpad=10
        dcolor=0.01
        pos1 = axx.get_position()
        cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
        cbar=plt.colorbar(mappable=c2,cax=cax,orientation="horizontal",ticklocation='top',aspect=35)
        cbar.set_label("SST [°C]",labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)
        axx.set_label('')
        axx.set_xlabel('')
        axx.set_ylabel('')
        axx.set_title('')
        
        axx.xaxis.set_ticks([130,140,150,160,170]) 
        labels = [item.get_text()+'°E' for item in axx.get_xticklabels()]
        axx.set_xticklabels(labels)
        
        axx.yaxis.set_ticks([30,35,40,45,50]) 
        labels = [item.get_text()+'°N' for item in axx.get_yticklabels()]
        axx.set_yticklabels(labels)
        
        #Typical size
        dx=1000
        Rt=6356.7523e3
        conv=(np.pi*Rt*1e-3)/180
        lat1,lat2=U10Ml.lat.min().data,U10Ml.lat.min().data+1*dx/conv
        lon1,lon2=U10Ml.lon.min().data,U10Ml.lon.min().data+1*dx/conv

        arr=axx.annotate('',size=15,xy=(lon1,lat1), xytext=(lon2,lat1),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='k', lw=6),annotation_clip=False)
        arr=axx.annotate('',size=15,xy=(lon1,lat1), xytext=(lon1,lat2),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='k', lw=6),annotation_clip=False)
        axx.text(lon1,28,str(dx)+' km',size=15,zorder=8,horizontalalignment='center',verticalalignment='bottom',\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        axx.text((lon2+lon1)/2,(lat1+0.4),str(dx)+' km',size=15,zorder=8,horizontalalignment='center',verticalalignment='bottom',\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
        
        # --------------------------------------------------------------
        axx=plt.subplot(gs[0,1])
        dw=0
        
        axx.spines['bottom'].set_color(co)
        axx.spines['bottom'].set_linewidth(linewidth+dw)
        axx.spines['top'].set_color(co)
        axx.spines['top'].set_linewidth(linewidth+dw)
        axx.spines['left'].set_color(co)
        axx.spines['left'].set_linewidth(linewidth+dw)
        axx.spines['right'].set_color(co)
        axx.spines['right'].set_linewidth(linewidth+dw)
        
        axx.text(tx,ty,'b',weight="bold",size=25,horizontalalignment='center',verticalalignment='center',\
                transform = axx.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
        loni=[148,157]
        lati=[38,41.8]
        co='r'
        linestyle='solid'
        linewidth=4
        axx.vlines(loni[0],ymin=lati[0],ymax=lati[1],colors=co,linestyles=linestyle,linewidths=linewidth,zorder=5)
        axx.vlines(loni[1],ymin=lati[0],ymax=lati[1],colors=co,linestyles=linestyle,linewidths=linewidth,zorder=5)
        axx.hlines(lati[0],xmin=loni[0],xmax=loni[1],colors=co,linestyles=linestyle,linewidths=linewidth,zorder=5)
        axx.hlines(lati[1],xmin=loni[0],xmax=loni[1],colors=co,linestyles=linestyle,linewidths=linewidth,zorder=5)
        co='b'
        linewidth=4
        
        sstlev=np.arange(0,24,0.75)
        
        c2=(EFLUX.EFLUX).plot.contourf(ax=axx,\
                levels=nlev,cmap='coolwarm',alpha=1,add_colorbar=False,extend='max')
                
        c3=(SST2.Theta).plot.contour(ax=axx,\
                levels=sstlev,linewidths=1.3,colors='k',alpha=1)

        lkey=15
        alpha=18 #14
        quiverscale=2.5e2
        quiverwidth=0.004 
        
        lonnn,lattt=np.meshgrid(U10M.lon.data,U10M.lat.data)     
        q=axx.quiver(lonnn[::alpha,::alpha],lattt[::alpha,::alpha],\
                U10M.U10M.data[::alpha,::alpha],V10M.V10M.data[::alpha,::alpha],zorder=2,\
                width=quiverwidth,scale=quiverscale)#,angles='xy', scale_units='xy')
        
        txq,tyq=0.86,0.96
        tt=axx.text(txq+0.033,tyq,'                 ',weight="bold",size=20,horizontalalignment='center',verticalalignment='center',\
                transform = axx.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        tt.set_zorder(2)
        qt=axx.quiverkey(q, X=txq, Y=tyq, U=lkey,
             label=str(lkey)+' m.s$^{-1}$', labelpos='E',fontproperties={'weight': 'bold','size': 15},zorder=3)#,bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        t = qt.text.set_zorder(3)


        pos1 = axx.get_position()
        cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
        cbar=plt.colorbar(mappable=c2,cax=cax,orientation="horizontal",ticklocation='top',aspect=35)
        cbar.set_label("LHF [W.m$^{-2}$]",labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)
        
        axx.set_xlabel('')
        axx.set_ylabel('')
        
        axx.xaxis.set_ticks([145,150,155,160]) 
        labels = [item.get_text()+'°E' for item in axx.get_xticklabels()]
        axx.set_xticklabels(labels)

        axx.yaxis.set_ticks([32,34,36,38,40]) 
        labels = [item.get_text()+'°N' for item in axx.get_yticklabels()]
        axx.set_yticklabels(labels)
        
        dx=400
        Rt=6356.7523e3
        conv=(np.pi*Rt*1e-3)/180
        lat1,lat2=U10M.lat.min().data,U10M.lat.min().data+1*dx/conv
        lon1,lon2=U10M.lon.min().data,U10M.lon.min().data+1*dx/conv

        arr=axx.annotate('',size=15,xy=(lon1,lat1), xytext=(lon2,lat1),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='k', lw=6),annotation_clip=False)
        arr=axx.annotate('',size=15,xy=(lon1,lat1), xytext=(lon1,lat2),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='k', lw=6),annotation_clip=False)
        axx.text(lon1,(lat2+lat1)/2,str(dx)+' km',size=15,zorder=8,horizontalalignment='center',verticalalignment='bottom',\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        axx.text((lon2+lon1)/2,(lat1+0.15),str(dx)+' km',size=15,zorder=8,horizontalalignment='center',verticalalignment='bottom',\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
        
        # --------------------------------------------------------------
        axx1= plt.subplot(gs[1,1],sharex=axx,sharey=axx)
        axx1.spines['bottom'].set_color(co)
        axx1.spines['bottom'].set_linewidth(linewidth+dw)
        axx1.spines['top'].set_color(co)
        axx1.spines['top'].set_linewidth(linewidth+dw)
        axx1.spines['left'].set_color(co)
        axx1.spines['left'].set_linewidth(linewidth+dw)
        axx1.spines['right'].set_color(co)
        axx1.spines['right'].set_linewidth(linewidth+dw)
        
        axx1.text(tx,ty,'d',weight="bold",size=25,horizontalalignment='center',verticalalignment='center',\
                transform = axx1.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
        
        levels=np.arange(-40,41,2)
        maxi=(divU10.divU10*1e5).max()
        c2=(divU10.divU10*1e5).plot.contourf(ax=axx1,\
                levels=levels,cmap='seismic',extend='both',add_colorbar=False)
        
        c3=(SST2.Theta).plot.contour(ax=axx1,\
                levels=sstlev,linewidths=1.3,colors='k',alpha=1)

        lkey=15
        
        lonnn,lattt=np.meshgrid(U10M.lon.data,U10M.lat.data)     
        q=axx1.quiver(lonnn[::alpha,::alpha],lattt[::alpha,::alpha],\
                U10M.U10M.data[::alpha,::alpha],V10M.V10M.data[::alpha,::alpha],zorder=2,\
                width=quiverwidth,scale=quiverscale)#,angles='xy', scale_units='xy')
        
        txq,tyq=0.86,0.96
        tt=axx1.text(txq+0.033,tyq,'                 ',weight="bold",size=20,horizontalalignment='center',verticalalignment='center',\
                transform = axx1.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        tt.set_zorder(2)
        qt=axx1.quiverkey(q, X=txq, Y=tyq, U=lkey,
             label=str(lkey)+' m.s$^{-1}$', labelpos='E',fontproperties={'weight': 'bold','size': 15},zorder=3)#,bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        t = qt.text.set_zorder(3)
        
        pos1 = axx1.get_position()
        cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
        cbar=plt.colorbar(mappable=c2,cax=cax,orientation="horizontal",ticklocation='top',aspect=35)
        cbar.set_label('$\mathbf{\\nabla}. \mathbf{u_{10}}$ [m.s$^{-1}$.(100km)$^{-1}$]',labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)
        
        axx1.set_ylabel('')
        axx1.set_xlabel('')
        
        dx=400
        Rt=6356.7523e3
        conv=(np.pi*Rt*1e-3)/180
        lat1,lat2=U10M.lat.min().data,U10M.lat.min().data+1*dx/conv
        lon1,lon2=U10M.lon.min().data,U10M.lon.min().data+1*dx/conv

        arr=axx1.annotate('',size=15,xy=(lon1,lat1), xytext=(lon2,lat1),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='k', lw=6),annotation_clip=False)
        arr=axx1.annotate('',size=15,xy=(lon1,lat1), xytext=(lon1,lat2),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='k', lw=6),annotation_clip=False)
        axx1.text(lon1,(lat2+lat1)/2,str(dx)+' km',size=15,zorder=8,horizontalalignment='center',verticalalignment='bottom',\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        axx1.text((lon2+lon1)/2,(lat1+0.15),str(dx)+' km',size=15,zorder=8,horizontalalignment='center',verticalalignment='bottom',\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
        # --------------------------------------------------------------
        axx1= plt.subplot(gs[1,0])
        axx1.spines['bottom'].set_color(co)
        axx1.spines['bottom'].set_linewidth(linewidth+dw)
        axx1.spines['top'].set_color(co)
        axx1.spines['top'].set_linewidth(linewidth+dw)
        axx1.spines['left'].set_color(co)
        axx1.spines['left'].set_linewidth(linewidth+dw)
        axx1.spines['right'].set_color(co)
        axx1.spines['right'].set_linewidth(linewidth+dw)
        
        axx1.text(tx,ty,'c',weight="bold",size=25,horizontalalignment='center',verticalalignment='center',\
                transform = axx1.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
                
        levels=np.arange(-0.75,0.75,0.025)

        LON,LAT=np.meshgrid(vort.lon.data,vort.lat.data)
        f=2*7.2921e-5*np.sin(LAT*np.pi/180)
        c2=(vort.vort/f).plot.contourf(ax=axx1,\
                levels=levels,cmap=cmocean.cm.balance,extend='both',add_colorbar=False)
        
        pos1 = axx1.get_position()
        cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
        cbar=plt.colorbar(mappable=c2,cax=cax,ticks=[-0.75,-0.5,-0.25,-0,0.25,0.5,0.75],orientation="horizontal",ticklocation='top',aspect=35)
        cbar.set_label('$\zeta /f$',labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)
        
        lkey=15
        
        lonnn,lattt=np.meshgrid(U10M.lon.data,U10M.lat.data)     
        q=axx1.quiver(lonnn[::alpha,::alpha],lattt[::alpha,::alpha],\
                U10M.U10M.data[::alpha,::alpha],V10M.V10M.data[::alpha,::alpha],zorder=2,\
                width=quiverwidth,scale=quiverscale)
        
        txq,tyq=0.86,0.96
        tt=axx1.text(txq+0.033,tyq,'                 ',weight="bold",size=20,horizontalalignment='center',verticalalignment='center',\
                transform = axx1.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        tt.set_zorder(2)
        qt=axx1.quiverkey(q, X=txq, Y=tyq, U=lkey,
             label=str(lkey)+' m.s$^{-1}$', labelpos='E',fontproperties={'weight': 'bold','size': 15},zorder=3)#,bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        t = qt.text.set_zorder(3)
        
        
        axx1.set_xlabel('')
        axx1.set_ylabel('')
        
        axx1.xaxis.set_ticks([145,150,155,160]) 
        labels = [item.get_text()+'°E' for item in axx1.get_xticklabels()]
        axx1.set_xticklabels(labels)

        axx1.yaxis.set_ticks([32,34,36,38,40]) 
        labels = [item.get_text()+'°N' for item in axx1.get_yticklabels()]
        axx1.set_yticklabels(labels)
        
        dx=400
        Rt=6356.7523e3
        conv=(np.pi*Rt*1e-3)/180
        lat1,lat2=U10M.lat.min().data,U10M.lat.min().data+1*dx/conv
        lon1,lon2=U10M.lon.min().data,U10M.lon.min().data+1*dx/conv

        arr=axx1.annotate('',size=15,xy=(lon1,lat1), xytext=(lon2,lat1),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='k', lw=6),annotation_clip=False)
        arr=axx1.annotate('',size=15,xy=(lon1,lat1), xytext=(lon1,lat2),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='k', lw=6),annotation_clip=False)
        axx1.text(lon1,(lat2+lat1)/2,str(dx)+' km',size=15,zorder=8,horizontalalignment='center',verticalalignment='bottom',\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        axx1.text((lon2+lon1)/2,(lat1+0.15),str(dx)+' km',size=15,zorder=8,horizontalalignment='center',verticalalignment='bottom',\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))   
        

# ----------------- Plot
date=datetime(2020,12,27,13)
endfile=datetime.strftime(date,'_2D_'+'%Y%m%d_%H%M' + '.nc')
plott(endfile)
plt.savefig('/path/to/your/plots/fig2.png')
plt.close()


