# Coded by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program to plot Extended data Fig.3

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec


def plott(endname,endname2):
        
        path='/path/to/your/extracted_data/2D/'
        
        ################## 5-day rms ###############################
        lonmin,lonmax=148,157
        latmin,latmax=38,42
        name='U10M'
        
        file=path+name+'_mean'+endname[4:]
        U10M=xr.open_dataset(file).sel(lon=slice(lonmin,lonmax)).sel(lat=slice(latmin,latmax))
        
        name='V10M'
        file=path+name+'_mean'+endname[4:]
        V10M=xr.open_dataset(file).sel(lon=slice(lonmin,lonmax)).sel(lat=slice(latmin,latmax))

        name='Theta'
        file=path+name+'_mean'+endname[4:]
        SST2=xr.open_dataset(file).sel(lon=slice(lonmin,lonmax)).sel(lat=slice(latmin,latmax))
        
        name='CNVMFC_3000m'
        file=path+name+endname
        CNVMFC_3000m=xr.open_dataset(file).sel(lon=slice(lonmin,lonmax)).sel(lat=slice(latmin,latmax))

        name='DTHDTCN_3000m'
        file=path+name+endname
        DTHDTCN_3000m=xr.open_dataset(file).sel(lon=slice(lonmin,lonmax)).sel(lat=slice(latmin,latmax))
        
        name='PRECCON'
        file=path+name+endname
        PRECCON=xr.open_dataset(file).sel(lon=slice(lonmin,lonmax)).sel(lat=slice(latmin,latmax))

        size=20
        parameters = {'axes.labelsize': size, 'axes.titlesize': size,"axes.labelsize" : size,"xtick.labelsize" : size,\
                "ytick.labelsize" : size,"ytick.major.size" : 5,"ytick.major.width" : 2,"ytick.minor.size" : 3,"ytick.minor.width" : 1,\
       "xtick.major.size" : 5,"xtick.major.width" : 2,"xtick.minor.size" : 3,"xtick.minor.width" : 1}
        plt.rcParams.update(parameters)

        fig=plt.figure(figsize=(25,17))
        
        gs=gridspec.GridSpec(2,3)#,height_ratios= [2,1,1])
        gs.update(bottom=0.03, top=0.92, left=0.04, right=0.98,
                        wspace=0.06, hspace=0.25)
        alphax=5
        alphay=9
        quiverscale=1.3e2 #1.5e2
        quiverwidth=0.004
        
        dcolor=0.008
        
        panelsize=22
        tx,ty=0.04,0.95
        labelpad=10
        sstlevels=np.arange(5,21,1)
        tick_font_size=18
        # --------------------------------------------------------------
        axx50=plt.subplot(gs[0,0])
        # plt.setp(axx50.get_yticklabels(), visible=False)
        # plt.setp(axx50.get_xticklabels(), visible=False)
        axx50.text(tx,ty,'a',weight="bold",size=panelsize,horizontalalignment='center',verticalalignment='center',\
                transform = axx50.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
        c3=(SST2.Theta).plot.contour(ax=axx50,\
                levels=sstlevels,linewidths=2,colors='k',alpha=0.8)
        axx50.clabel(c3, inline=1, fontsize=20)
        
        max=2000
        # norm = colors.Normalize(vmin=-max, vmax=max)
        maxi=(CNVMFC_3000m.CNVMFC*3600*24).max()
        print(maxi)
        c2=(CNVMFC_3000m.CNVMFC*3600*24).plot.contourf(ax=axx50,\
                levels=np.arange(0,6001,200),cmap='YlOrRd',extend='max',add_colorbar=False)#,norm=norm)
                
        pos1 = axx50.get_position()
        cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
        cbar=plt.colorbar(mappable=c2,cax=cax,orientation="horizontal",ticklocation='top',aspect=35)
        cbar.set_label("Convective Mass Flux at $3$ km [kg.m$^{-2}$.day$^{-1}$]",labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)
        
        
        lonnn,lattt=np.meshgrid(U10M.lon.data,U10M.lat.data)
        
        q=axx50.quiver(lonnn[::alphax,::alphay],lattt[::alphax,::alphay],\
                U10M.U10M.data[::alphax,::alphay],V10M.V10M.data[::alphax,::alphay],zorder=2,\
                width=quiverwidth,scale=quiverscale)
        
        lkey=5
        txq,tyq=0.84,0.96

        tt=axx50.text(txq+0.05,tyq,'              ',weight="bold",size=22,horizontalalignment='center',verticalalignment='center',\
                transform = axx50.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        tt.set_zorder(2)
        qt=axx50.quiverkey(q, X=txq, Y=tyq, U=lkey,
             label=str(lkey)+' m.s$^{-1}$', labelpos='E',fontproperties={'weight': 'bold','size': 15},zorder=3)#,bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        t = qt.text.set_zorder(3)
        
        axx50.set_xlabel('')
        axx50.set_ylabel('')
        
        axx50.yaxis.set_ticks([38,39,40,41]) 
        labels = [item.get_text()+'°N' for item in axx50.get_yticklabels()]
        axx50.set_yticklabels(labels)
         
        axx50.xaxis.set_ticks([149,151,153,155,157])  #150,152,154,156
        labels = [item.get_text()+'°E' for item in axx50.get_xticklabels()]
        axx50.set_xticklabels(labels)  
        
        dx=200
        Rt=6356.7523e3
        conv=(np.pi*Rt*1e-3)/180
        lat1,lat2=U10M.lat.min().data,U10M.lat.min().data+1*dx/conv
        lon1,lon2=U10M.lon.min().data,U10M.lon.min().data+1*dx/conv
        # ypos=60*hscale
        arr=axx50.annotate('',size=15,xy=(lon1,lat1), xytext=(lon2,lat1),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='k', lw=6),annotation_clip=False)
        arr=axx50.annotate('',size=15,xy=(lon1,lat1), xytext=(lon1,lat2),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='k', lw=6),annotation_clip=False)
        axx50.text(lon1,(lat2+lat1)/2-0.3,str(dx)+' km',size=15,zorder=8,horizontalalignment='center',verticalalignment='bottom',\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        axx50.text((lon2+lon1)/2,(lat1+0.08),str(dx)+' km',size=15,zorder=8,horizontalalignment='center',verticalalignment='bottom',\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
        # --------------------------------------------------------------
        axx50=plt.subplot(gs[0,1],sharex=axx50,sharey=axx50)
        plt.setp(axx50.get_yticklabels(), visible=False)
        # plt.setp(axx50.get_xticklabels(), visible=False)
        axx50.text(tx,ty,'b',weight="bold",size=panelsize,horizontalalignment='center',verticalalignment='center',\
                transform = axx50.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))


        c3=(SST2.Theta).plot.contour(ax=axx50,\
                levels=sstlevels,linewidths=2,colors='k',alpha=0.8)
        axx50.clabel(c3, inline=1, fontsize=20)
        
        c2=(DTHDTCN_3000m.DTHDTCN*3600*24).plot.contourf(ax=axx50,\
        levels=np.arange(0,30,2),cmap='YlOrRd',extend='max',add_colorbar=False)
        
        pos1 = axx50.get_position()
        cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
        cbar=plt.colorbar(mappable=c2,cax=cax,orientation="horizontal",ticklocation='top',aspect=35)
        cbar.set_label("Convective diabatic heating at $3$ km [K.day$^{-1}$]",labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)
        
        lonnn,lattt=np.meshgrid(U10M.lon.data,U10M.lat.data)
        
        axx50.quiver(lonnn[::alphax,::alphay],lattt[::alphax,::alphay],\
                U10M.U10M.data[::alphax,::alphay],V10M.V10M.data[::alphax,::alphay],zorder=2,\
                width=quiverwidth,scale=quiverscale)
        axx50.set_ylabel('')
        axx50.set_xlabel('')
        
        # --------------------------------------------------------------
        axx50=plt.subplot(gs[0,2],sharex=axx50,sharey=axx50)
        plt.setp(axx50.get_yticklabels(), visible=False)
        axx50.text(tx,ty,'c',weight="bold",size=panelsize,horizontalalignment='center',verticalalignment='center',\
                transform = axx50.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        

        c3=(SST2.Theta).plot.contour(ax=axx50,\
                levels=sstlevels,linewidths=2,colors='k',alpha=0.8)
        axx50.clabel(c3, inline=1, fontsize=20)
        
        levels=np.arange(0,41,1)
        maxi=(PRECCON.PRECCON*3600*24).max()
        print(maxi)
        c2=(PRECCON.PRECCON*3600*24).plot.contourf(ax=axx50,\
                levels=levels,cmap='Blues',extend='max',add_colorbar=False)     
                

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
        
        ################## 2-month rms ###############################
        endname=endname2
        
        name='U10M'
        file=path+name+'_mean'+endname[4:]
        U10M=xr.open_dataset(file)#.sel(lon=slice(lonmin,lonmax)).sel(lat=slice(latmin,latmax))
        
        name='V10M'
        file=path+name+'_mean'+endname[4:]
        V10M=xr.open_dataset(file)#.sel(lon=slice(lonmin,lonmax)).sel(lat=slice(latmin,latmax))

        
        name='CNVMFC_3000m'
        file=path+name+endname
        CNVMFC_3000m=xr.open_dataset(file)
        
        name='DTHDTCN_3000m'
        file=path+name+endname
        DTHDTCN_3000m=xr.open_dataset(file)

        name='PRECCON'
        file=path+name+endname
        PRECCON=xr.open_dataset(file)

        name='Theta'
        file=path+name+'_mean'+endname[4:]
        SST2=xr.open_dataset(file)#.sel(lon=slice(lonmin,lonmax)).sel(lat=slice(latmin,latmax))
        
        alphax,alphay=11,17
        quiverscale=1.7e2 #1.5e2
        quiverwidth=0.0035
        
        dcolor=0.008
        
        panelsize=22
        tx,ty=0.04,0.95
        labelpad=10
        lkey=10
        txq,tyq=0.83,0.96
        
        # --------------------------------------------------------------
        axx50=plt.subplot(gs[1,0])
        # plt.setp(axx50.get_yticklabels(), visible=False)
        # plt.setp(axx50.get_xticklabels(), visible=False)
        axx50.text(tx,ty,'d',weight="bold",size=panelsize,horizontalalignment='center',verticalalignment='center',\
                transform = axx50.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
        c3=(SST2.Theta).plot.contour(ax=axx50,\
                levels=sstlevels,linewidths=2,colors='k',alpha=0.8)
        axx50.clabel(c3, inline=1, fontsize=20)
        
        # norm = colors.Normalize(vmin=-1000, vmax=1000)
        # maxi=(CNVMFC_3000m.CNVMFC*3600*24).max()
        levels=30
        c2=(CNVMFC_3000m.CNVMFC*3600*24).plot.contourf(ax=axx50,\
                levels=levels,cmap='YlOrRd',extend='max',add_colorbar=False)#,norm=norm)
                
        pos1 = axx50.get_position()
        cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
        cbar=plt.colorbar(mappable=c2,cax=cax,orientation="horizontal",ticklocation='top',aspect=35)
        cbar.set_label("Convective Mass Flux at 3 km [kg.m$^{-2}$.day$^{-1}$]",labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)
        
        
        lonnn,lattt=np.meshgrid(U10M.lon.data,U10M.lat.data)
        
        q=axx50.quiver(lonnn[::alphax,::alphay],lattt[::alphax,::alphay],\
                U10M.U10M.data[::alphax,::alphay],V10M.V10M.data[::alphax,::alphay],zorder=2,\
                width=quiverwidth,scale=quiverscale)
        
        lkey=10
        txq,tyq=0.83,0.96
        tt=axx50.text(txq+0.045,tyq,'                ',weight="bold",size=22,horizontalalignment='center',verticalalignment='center',\
                transform = axx50.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        tt.set_zorder(2)
        qt=axx50.quiverkey(q, X=txq, Y=tyq, U=lkey,transform = axx50.transAxes,
             label=str(lkey)+' m.s$^{-1}$', labelpos='E',fontproperties={'weight': 'bold','size': 15},zorder=3)#,bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        t = qt.text.set_zorder(3)
        
        axx50.set_xlabel('')
        axx50.set_ylabel('')
        
        axx50.yaxis.set_ticks([32,34,36,38,40]) 
        labels = [item.get_text()+'°N' for item in axx50.get_yticklabels()]
        axx50.set_yticklabels(labels)
        
        axx50.xaxis.set_ticks([145,150,155,160]) 
        labels = [item.get_text()+'°E' for item in axx50.get_xticklabels()]
        axx50.set_xticklabels(labels)  
        
        dx=400
        Rt=6356.7523e3
        conv=(np.pi*Rt*1e-3)/180
        lat1,lat2=U10M.lat.min().data,U10M.lat.min().data+1*dx/conv
        lon1,lon2=U10M.lon.min().data,U10M.lon.min().data+1*dx/conv
        # ypos=60*hscale
        arr=axx50.annotate('',size=15,xy=(lon1,lat1), xytext=(lon2,lat1),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='k', lw=6),annotation_clip=False)
        arr=axx50.annotate('',size=15,xy=(lon1,lat1), xytext=(lon1,lat2),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='k', lw=6),annotation_clip=False)
        axx50.text(lon1,(lat2+lat1)/2,str(dx)+' km',size=15,zorder=8,horizontalalignment='center',verticalalignment='bottom',\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        axx50.text((lon2+lon1)/2,(lat1+0.15),str(dx)+' km',size=15,zorder=8,horizontalalignment='center',verticalalignment='bottom',\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
        # --------------------------------------------------------------
        axx50=plt.subplot(gs[1,1],sharex=axx50)
        plt.setp(axx50.get_yticklabels(), visible=False)
        # plt.setp(axx50.get_xticklabels(), visible=False)
        axx50.text(tx,ty,'e',weight="bold",size=panelsize,horizontalalignment='center',verticalalignment='center',\
                transform = axx50.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))


        c3=(SST2.Theta).plot.contour(ax=axx50,\
                levels=sstlevels,linewidths=2,colors='k',alpha=0.8)
        axx50.clabel(c3, inline=1, fontsize=20)

        levels=np.arange(0,20,1)
        c2=(DTHDTCN_3000m.DTHDTCN*3600*24).plot.contourf(ax=axx50,\
                levels=levels,cmap='YlOrRd',extend='max',add_colorbar=False)
        
        pos1 = axx50.get_position()
        cax = fig.add_axes([pos1.x0, pos1.y1+dcolor, pos1.x1-pos1.x0, dcolor])
        cbar=plt.colorbar(mappable=c2,cax=cax,orientation="horizontal",ticklocation='top',aspect=35)
        cbar.set_label("Convective diabatic heating at $3$ km [K.day$^{-1}$]",labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)
        
        
        lonnn,lattt=np.meshgrid(U10M.lon.data,U10M.lat.data)
        
        axx50.quiver(lonnn[::alphax,::alphay],lattt[::alphax,::alphay],\
                U10M.U10M.data[::alphax,::alphay],V10M.V10M.data[::alphax,::alphay],zorder=2,\
                width=quiverwidth,scale=quiverscale)
        axx50.set_ylabel('')
        axx50.set_xlabel('')
        
        
        # --------------------------------------------------------------
        axx50=plt.subplot(gs[1,2],sharex=axx50)
        plt.setp(axx50.get_yticklabels(), visible=False)
        axx50.text(tx,ty,'f',weight="bold",size=panelsize,horizontalalignment='center',verticalalignment='center',\
                transform = axx50.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        

        c3=(SST2.Theta).plot.contour(ax=axx50,\
                levels=sstlevels,linewidths=2,colors='k',alpha=0.8)
        axx50.clabel(c3, inline=1, fontsize=20)
        
        levels=np.arange(0,20,1)
        maxi=(PRECCON.PRECCON*3600*24).max()
        c2=(PRECCON.PRECCON*3600*24).plot.contourf(ax=axx50,\
                levels=levels,cmap='Blues',extend='max',add_colorbar=False)     
                

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


plott('_rms_16dec21dec_NW0sup_2D.nc','_rms_decjan_NW0sup_2D.nc')
plt.savefig('/path/to/your/plots/ExtDataFig3.png')

plt.close()


