# Coded by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program to plot Extended data Fig.7

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

def plott(endname):
        
        path='/path/to/your/extracted_data/2D/'
        name='U10M'
        file=path+name+endname
        U10M=xr.open_dataset(file)
        
        name='V10M'
        file=path+name+endname
        V10M=xr.open_dataset(file)

        name='PBLH'
        file=path+name+endname
        PBLH=xr.open_dataset(file)
        
        name='PBLH'
        file=path+name+endname
        PBLH=xr.open_dataset(file)
        
        name='MFC'
        file=path+name+endname
        MFC=xr.open_dataset(file)
        
        name='CNVMFC_2000m'
        file=path+name+endname
        CNVMFC_2000m=xr.open_dataset(file)
        
        name='DTHDTCN_2000m'
        file=path+name+endname
        DTHDTCN_2000m=xr.open_dataset(file)
        
        name='PRECCON'
        file=path++name+endname
        PRECCON=xr.open_dataset(file)

        name='gradTheta_dwind'
        file=path+name+endname
        gradSST_dwind2=xr.open_dataset(file)
        
        name='CAPE'
        file=path+name+endname
        CAPE=xr.open_dataset(file)
        
        name='EFLUX'
        file=path+name+endname
        EFLUX=xr.open_dataset(file)
        
        name='Theta'
        file=path+name+endname
        SST2=xr.open_dataset(file)
    
        name='divU10'
        file=path+name+endname
        divU10=xr.open_dataset(file)
        
        size=20
        parameters = {'axes.labelsize': size, 'axes.titlesize': size,"axes.labelsize" : size,"xtick.labelsize" : size,\
                "ytick.labelsize" : size,"ytick.major.size" : 5,"ytick.major.width" : 2,"ytick.minor.size" : 3,"ytick.minor.width" : 1,\
       "xtick.major.size" : 5,"xtick.major.width" : 2,"xtick.minor.size" : 3,"xtick.minor.width" : 1}
        plt.rcParams.update(parameters)

        fig=plt.figure(figsize=(25,25))
        
        gs=gridspec.GridSpec(3,3)
        gs.update(bottom=0.03, top=0.95, left=0.04, right=0.98,
                        wspace=0.06, hspace=0.22)

        alphax,alphay=11,17
        sstlev=10
        quiverscale=1.7e2 #1.5e2
        quiverwidth=0.0035
        
        dcolor=0.008
        
        panelsize=22
        tx,ty=0.04,0.95
        labelpad=10
        
        # --------------------------------------------------------------
        axx=plt.subplot(gs[0,0])
        plt.setp(axx.get_xticklabels(), visible=False)
        # plt.setp(axx.get_yticklabels(), visible=False)
        
        axx.text(tx,ty,'a',weight="bold",size=panelsize,horizontalalignment='center',verticalalignment='center',\
                transform = axx.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
        levels=np.arange(50,360,10)
        c2=(EFLUX.EFLUX).plot.contourf(ax=axx,\
                levels=levels,cmap='coolwarm',alpha=1,add_colorbar=False,extend='max')
                
        sstlevels=np.arange(5,21,1)#[20,17.5,15,12.5,10,7.5,5]
        c3=(SST2.Theta).plot.contour(ax=axx,\
                levels=sstlevels,linewidths=2,colors='k',alpha=0.8)
        axx.clabel(c3, inline=1, fontsize=20)

        lonnn,lattt=np.meshgrid(U10M.lon.data,U10M.lat.data)     
        q=axx.quiver(lonnn[::alphax,::alphay],lattt[::alphax,::alphay],\
                U10M.U10M.data[::alphax,::alphay],V10M.V10M.data[::alphax,::alphay],zorder=2,\
                width=quiverwidth,scale=quiverscale)
        
        lkey=10
        txq,tyq=0.83,0.96
        tt=axx.text(txq+0.045,tyq,'                ',weight="bold",size=22,horizontalalignment='center',verticalalignment='center',\
                transform = axx.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        tt.set_zorder(2)
        qt=axx.quiverkey(q, X=txq, Y=tyq, U=lkey,
             label=str(lkey)+' m.s$^{-1}$', labelpos='E',fontproperties={'weight': 'bold','size': 15},zorder=3)#,bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        t = qt.text.set_zorder(3)

        tick_font_size=18
        pos1 = axx.get_position()
        cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
        cbar=plt.colorbar(mappable=c2,cax=cax,orientation="horizontal",ticklocation='top',aspect=35)
        cbar.set_label("LHF [W.m$^{-2}$]",labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)
        
        axx.set_xlabel('')
        axx.set_ylabel('')
        
        axx.yaxis.set_ticks([32,34,36,38,40])
        axx.xaxis.set_ticks([145,150,155,160]) 
        labels = [item.get_text()+'°N' for item in axx.get_yticklabels()]
        axx.set_yticklabels(labels)
        
        
        lonsel=152.5
        ymin,ymax=35,40.5
        axx.vlines(lonsel,ymin=ymin,ymax=ymax,colors='b',linestyles='dashed',linewidths=5)
        lonsel=152.5
        ymin,ymax=40,41.5
        axx.vlines(lonsel,ymin=ymin,ymax=ymax,colors='r',linestyles='dashed',linewidths=5)
        
        loni=[148,157]
        lati=[37.5,41.8]
        co='k'
        linestyle='solid'
        linewidth=4
        axx.vlines(loni[0],ymin=lati[0],ymax=lati[1],colors=co,linestyles=linestyle,linewidths=linewidth)
        axx.vlines(loni[1],ymin=lati[0],ymax=lati[1],colors=co,linestyles=linestyle,linewidths=linewidth)
        axx.hlines(lati[0],xmin=loni[0],xmax=loni[1],colors=co,linestyles=linestyle,linewidths=linewidth)
        axx.hlines(lati[1],xmin=loni[0],xmax=loni[1],colors=co,linestyles=linestyle,linewidths=linewidth)
        
        axx.scatter(lonsel,36.5,marker='*',color='b',s=1000,clip_on=False,zorder=4)
        axx.scatter(lonsel,39.1,marker='*',color='b',s=1000,clip_on=False,zorder=4)        
        axx.scatter(lonsel,40.65,marker='*',color='r',s=1000,clip_on=False,zorder=4)

        dx=400
        Rt=6356.7523e3
        conv=(np.pi*Rt*1e-3)/180
        lat1,lat2=U10M.lat.min().data,U10M.lat.min().data+1*dx/conv
        lon1,lon2=U10M.lon.min().data,U10M.lon.min().data+1*dx/conv
        # ypos=60*hscale
        arr=axx.annotate('',size=15,xy=(lon1,lat1), xytext=(lon2,lat1),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='k', lw=6),annotation_clip=False)
        arr=axx.annotate('',size=15,xy=(lon1,lat1), xytext=(lon1,lat2),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='k', lw=6),annotation_clip=False)
        axx.text(lon1,(lat2+lat1)/2,str(dx)+' km',size=15,zorder=8,horizontalalignment='center',verticalalignment='bottom',\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        axx.text((lon2+lon1)/2,(lat1+0.15),str(dx)+' km',size=15,zorder=8,horizontalalignment='center',verticalalignment='bottom',\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        # --------------------------------------------------------------
        axx50=plt.subplot(gs[1,2],sharex=axx,sharey=axx)
        plt.setp(axx50.get_yticklabels(), visible=False)
        plt.setp(axx50.get_xticklabels(), visible=False)
        axx50.text(tx,ty,'f',weight="bold",size=panelsize,horizontalalignment='center',verticalalignment='center',\
                transform = axx50.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        

        c3=(SST2.Theta).plot.contour(ax=axx50,\
                levels=sstlevels,linewidths=2,colors='k',alpha=0.8)
        axx50.clabel(c3, inline=1, fontsize=20)
        
        # maxi=(W_2000mean.W*1e2).max()
        levels=np.arange(-25,25.1,1)
        c2=(MFC.MFC*1e3*3600*24).plot.contourf(ax=axx50,\
                levels=levels,cmap='seismic',extend='both',add_colorbar=False)
                
        pos1 = axx50.get_position()
        cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
        cbar=plt.colorbar(mappable=c2,cax=cax,orientation="horizontal",ticklocation='top',aspect=35)
        cbar.set_label("Moisture Flux Convergence at 10 m [g.kg$^{-1}$.day$^{-1}$]",labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)
        
        
        lonnn,lattt=np.meshgrid(U10M.lon.data,U10M.lat.data)
        
        axx50.quiver(lonnn[::alphax,::alphay],lattt[::alphax,::alphay],\
                U10M.U10M.data[::alphax,::alphay],V10M.V10M.data[::alphax,::alphay],zorder=2,\
                width=quiverwidth,scale=quiverscale)
        axx50.set_xlabel('')
        axx50.set_ylabel('')
        
        # --------------------------------------------------------------
        axx=plt.subplot(gs[0,1],sharex=axx,sharey=axx)
        plt.setp(axx.get_xticklabels(), visible=False)
        plt.setp(axx.get_yticklabels(), visible=False)
        
        axx.text(tx,ty,'b',weight="bold",size=panelsize,horizontalalignment='center',verticalalignment='center',\
                transform = axx.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))


        levels=np.arange(0,111,5)
        c2=(CAPE.CAPE).plot.contourf(ax=axx,\
                levels=levels,cmap='coolwarm',extend='max',add_colorbar=False)
        

        c3=(SST2.Theta).plot.contour(ax=axx,\
                levels=sstlevels,linewidths=2,colors='k',alpha=0.8)
        axx.clabel(c3, inline=1, fontsize=20)

        # print(((U10M.U10M**2 +V10M.V10M**2)**0.5).max())
        lonnn,lattt=np.meshgrid(U10M.lon.data,U10M.lat.data)     
        axx.quiver(lonnn[::alphax,::alphay],lattt[::alphax,::alphay],\
                U10M.U10M.data[::alphax,::alphay],V10M.V10M.data[::alphax,::alphay],zorder=2,\
                width=quiverwidth,scale=quiverscale)


        pos1 = axx.get_position()
        cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
        cbar=plt.colorbar(mappable=c2,cax=cax,orientation="horizontal",ticklocation='top',aspect=35)
        cbar.set_label("CAPE [J.kg$^{-1}$]",labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)

        axx.set_ylabel('')
        axx.set_xlabel('')
        
        
        # --------------------------------------------------------------
        axx=plt.subplot(gs[0,2])
        plt.setp(axx.get_xticklabels(), visible=False)
        plt.setp(axx.get_yticklabels(), visible=False)
        
        axx.text(tx,ty,'c',weight="bold",size=panelsize,horizontalalignment='center',verticalalignment='center',\
                transform = axx.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))


        levels=np.arange(1100,1801,25)
        c2=(PBLH.PBLH).plot.contourf(ax=axx,\
                levels=levels,cmap='coolwarm',extend='max',add_colorbar=False)
        
        c3=(SST2.Theta).plot.contour(ax=axx,\
                levels=sstlevels,linewidths=2,colors='k',alpha=0.8)
        axx.clabel(c3, inline=1, fontsize=20)

        # print(((U10M.U10M**2 +V10M.V10M**2)**0.5).max())
        lonnn,lattt=np.meshgrid(U10M.lon.data,U10M.lat.data)     
        axx.quiver(lonnn[::alphax,::alphay],lattt[::alphax,::alphay],\
                U10M.U10M.data[::alphax,::alphay],V10M.V10M.data[::alphax,::alphay],zorder=2,\
                width=quiverwidth,scale=quiverscale)


        pos1 = axx.get_position()
        cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
        cbar=plt.colorbar(mappable=c2,cax=cax,orientation="horizontal",ticklocation='top',aspect=35)
        cbar.set_label("PBLH [m]",labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)

        axx.set_ylabel('')
        axx.set_xlabel('')
        
        axx.yaxis.set_ticks([32,34,36,38,40]) 
        labels = [item.get_text()+'°N' for item in axx.get_yticklabels()]
        axx.set_yticklabels(labels)
        
        axx.xaxis.set_ticks([145,150,155,160]) 
        labels = [item.get_text()+'°E' for item in axx.get_xticklabels()]
        axx.set_xticklabels(labels)  
        
        # --------------------------------------------------------------
        axx5= plt.subplot(gs[1,0],sharex=axx) #plt.subplot(gs[1:,1])
        # plt.setp(axx5.get_yticklabels(), visible=False)
        plt.setp(axx5.get_xticklabels(), visible=False)
        axx5.text(tx,ty,'d',weight="bold",size=panelsize,horizontalalignment='center',verticalalignment='center',\
                transform = axx5.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
        levels=np.arange(-10,10.1,0.2)
        # maxi=(gradSST_dwind2.gradTheta_dwind*1e5).max()*0.9
        c2=(gradSST_dwind2.gradTheta_dwind*1e5).plot.contourf(ax=axx5,\
                levels=levels,cmap='seismic',extend='both',add_colorbar=False)
                

        c3=(SST2.Theta).plot.contour(ax=axx5,\
                levels=sstlevels,linewidths=2,colors='k',alpha=0.8)
        axx5.clabel(c3, inline=1, fontsize=20)
        
        lonnn,lattt=np.meshgrid(U10M.lon.data,U10M.lat.data)
        
        axx5.quiver(lonnn[::alphax,::alphay],lattt[::alphax,::alphay],\
                U10M.U10M.data[::alphax,::alphay],V10M.V10M.data[::alphax,::alphay],zorder=2,\
                width=quiverwidth,scale=quiverscale)

        axx5.set_ylabel('')
        axx5.set_xlabel('')
        
        axx5.xaxis.set_ticks([145,150,155,160]) 
        labels = [item.get_text()+'°N' for item in axx5.get_yticklabels()]
        axx5.set_yticklabels(labels)

        pos1 = axx5.get_position()
        cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
        cbar=plt.colorbar(mappable=c2,cax=cax,orientation="horizontal",ticklocation='top',aspect=35)
        cbar.set_label("$\mathbf{k}.\mathbf{\\nabla}$SST [°C.(100km)$^{-1}$]",labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)
        axx5.set_ylabel('')
        axx5.set_xlabel('')
        
        # --------------------------------------------------------------
        axx1= plt.subplot(gs[1,1],sharex=axx)
        plt.setp(axx1.get_xticklabels(), visible=False)
        plt.setp(axx1.get_yticklabels(), visible=False)
        axx1.text(tx,ty,'e',weight="bold",size=panelsize,horizontalalignment='center',verticalalignment='center',\
                transform = axx1.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
        levels=np.arange(-4,4.1,0.2)
        c2=(divU10.divU10*1e5).plot.contourf(ax=axx1,\
                levels=levels,cmap='seismic',extend='both',add_colorbar=False)
        

        c3=(SST2.Theta).plot.contour(ax=axx1,\
                levels=sstlevels,linewidths=2,colors='k',alpha=0.8)
        axx1.clabel(c3, inline=1, fontsize=20)

        lonnn,lattt=np.meshgrid(U10M.lon.data,U10M.lat.data)     
        axx1.quiver(lonnn[::alphax,::alphay],lattt[::alphax,::alphay],\
                U10M.U10M.data[::alphax,::alphay],V10M.V10M.data[::alphax,::alphay],zorder=2,\
                width=quiverwidth,scale=quiverscale)
        
        
        pos1 = axx1.get_position()
        cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
        cbar=plt.colorbar(mappable=c2,cax=cax,orientation="horizontal",ticklocation='top',aspect=35)
        cbar.set_label('$\mathbf{\\nabla}.\mathbf{u_{10}}$ [m.s$^{-1}$.(100km)$^{-1}$]',labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)
        
        axx1.set_ylabel('')
        axx1.set_xlabel('')
        
        # --------------------------------------------------------------
        axx50=plt.subplot(gs[2,0])
        # plt.setp(axx50.get_yticklabels(), visible=False)
        # plt.setp(axx50.get_xticklabels(), visible=False)
        axx50.text(tx,ty,'g',weight="bold",size=panelsize,horizontalalignment='center',verticalalignment='center',\
                transform = axx50.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
        c3=(SST2.Theta).plot.contour(ax=axx50,\
                levels=sstlevels,linewidths=2,colors='k',alpha=0.8)
        axx50.clabel(c3, inline=1, fontsize=20)
        
        # norm = colors.Normalize(vmin=-900, vmax=900)
        # maxi=(CNVMFC_3000m.CNVMFC*3600*24).max()
        # c2=(CNVMFC_3000m.CNVMFC*3600*24).plot.contourf(ax=axx50,\
        #         levels=np.arange(0,901,25),cmap='YlOrRd',extend='max',add_colorbar=False)#,norm=norm)
        
        c2=(CNVMFC_2000m.CNVMFC*3600*24).plot.contourf(ax=axx50,\
                levels=np.arange(0,1501,50),cmap='YlOrRd',extend='max',add_colorbar=False)#,norm=norm)
                
        pos1 = axx50.get_position()
        cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
        cbar=plt.colorbar(mappable=c2,cax=cax,orientation="horizontal",ticklocation='top',aspect=35)
        cbar.set_label("Convective Mass Flux at $2$ km [kg.m$^{-2}$.day$^{-1}$]",labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)
        
        
        lonnn,lattt=np.meshgrid(U10M.lon.data,U10M.lat.data)
        
        axx50.quiver(lonnn[::alphax,::alphay],lattt[::alphax,::alphay],\
                U10M.U10M.data[::alphax,::alphay],V10M.V10M.data[::alphax,::alphay],zorder=2,\
                width=quiverwidth,scale=quiverscale)
        axx50.set_xlabel('')
        axx50.set_ylabel('')
        
        
        axx50.yaxis.set_ticks([32,34,36,38,40]) 
        labels = [item.get_text()+'°N' for item in axx50.get_yticklabels()]
        axx50.set_yticklabels(labels)
        
        axx50.xaxis.set_ticks([145,150,155,160]) 
        labels = [item.get_text()+'°E' for item in axx50.get_xticklabels()]
        axx50.set_xticklabels(labels)  
        

        # --------------------------------------------------------------
        axx50=plt.subplot(gs[2,1],sharex=axx50)
        plt.setp(axx50.get_yticklabels(), visible=False)
        # plt.setp(axx50.get_xticklabels(), visible=False)
        axx50.text(tx,ty,'h',weight="bold",size=panelsize,horizontalalignment='center',verticalalignment='center',\
                transform = axx50.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))


        c3=(SST2.Theta).plot.contour(ax=axx50,\
                levels=sstlevels,linewidths=2,colors='k',alpha=0.8)
        axx50.clabel(c3, inline=1, fontsize=20)


        # c2=(DTHDTCN_3000m.DTHDTCN*3600*24).plot.contourf(ax=axx50,\
        #         levels=np.arange(0,3.1,0.25),cmap='YlOrRd',extend='max',add_colorbar=False)
        
        c2=(DTHDTCN_2000m.DTHDTCN*3600*24).plot.contourf(ax=axx50,\
                levels=np.arange(0,7.1,0.5),cmap='YlOrRd',extend='max',add_colorbar=False)
        
        pos1 = axx50.get_position()
        cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
        cbar=plt.colorbar(mappable=c2,cax=cax,orientation="horizontal",ticklocation='top',aspect=35)
        cbar.set_label("Convective diabatic heating at $2$ km [K.day$^{-1}$]",labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)
        
        lonnn,lattt=np.meshgrid(U10M.lon.data,U10M.lat.data)
        
        axx50.quiver(lonnn[::alphax,::alphay],lattt[::alphax,::alphay],\
                U10M.U10M.data[::alphax,::alphay],V10M.V10M.data[::alphax,::alphay],zorder=2,\
                width=quiverwidth,scale=quiverscale)
        axx50.set_ylabel('')
        axx50.set_xlabel('')
        
        # --------------------------------------------------------------
        axx50=plt.subplot(gs[2,2],sharex=axx50)
        plt.setp(axx50.get_yticklabels(), visible=False)
        axx50.text(tx,ty,'i',weight="bold",size=panelsize,horizontalalignment='center',verticalalignment='center',\
                transform = axx50.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        

        c3=(SST2.Theta).plot.contour(ax=axx50,\
                levels=sstlevels,linewidths=2,colors='k',alpha=0.8)
        axx50.clabel(c3, inline=1, fontsize=20)
        
        maxi=(PRECCON.PRECCON*3600*24).max()
        c2=(PRECCON.PRECCON*3600*24).plot.contourf(ax=axx50,\
                levels=np.arange(0,1.6,0.1),cmap='Blues',extend='max',add_colorbar=False)     
                

        pos1 = axx50.get_position()
        cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
        cbar=plt.colorbar(mappable=c2,cax=cax,orientation="horizontal",ticklocation='top',aspect=35)
        cbar.set_label("Convective precipitation [mm.day$^{-1}$]",labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)
        
        lonnn,lattt=np.meshgrid(U10M.lon.data,U10M.lat.data)
        
        axx50.quiver(lonnn[::alphax,::alphay],lattt[::alphax,::alphay],\
                U10M.U10M.data[::alphax,::alphay],V10M.V10M.data[::alphax,::alphay],zorder=2,\
                width=quiverwidth,scale=quiverscale)
        axx50.set_xlabel('')
        axx50.set_ylabel('')

plott('_mean_decjan_SW0_2D.nc')
plt.savefig('/path/to/your/plots/ExtDataFig7.png')

plt.close()


