# Coded by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program to plot Extended data Fig.6

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import cmocean

def plott(title,folder,lonsel,collection,collection2):
        zmax=4000
        rollat=6
        rolz=1
        latmin,latmax=35,40.5
        Rt=6356.7523e3
        
        #################### Warm sector ##################################
        name='H'
        H=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))
        
        name='wrho'
        wrho=xr.open_dataset(folder+name+collection).rolling(lat=rollat,depth=rolz,center=True).mean().sel(lat=slice(latmin,latmax))
        
        name='DTHDT_DELP'
        DTHDT_DELP=xr.open_dataset(folder+name+collection).rolling(lat=rollat,depth=rolz,center=True).mean().sel(lat=slice(latmin,latmax))
        
        name='P'
        P=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))
        
        name='PRECTOT'
        PRECTOT=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))
        
        name='CNVMFC'
        CNVMFC=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))
        
        name='QV'
        QV=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))
        
        name='W'
        W=xr.open_dataset(folder+name+collection).rolling(lat=rollat,depth=rolz,center=True).mean().sel(lat=slice(latmin,latmax))

        name='V'
        V=xr.open_dataset(folder+name+collection).rolling(lat=rollat,depth=rolz,center=True).mean().sel(lat=slice(latmin,latmax))

        divV=V.differentiate("lat")*180/(np.pi*Rt)

        name='Theta_ocean'
        SST=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))
       
        name='EFLUX'
        EFLUX=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))
        
        name='CAPE'
        CAPE=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))

        name='W_2000mean'
        W_2000mean=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))
        
        name='PBLH'
        PBLH=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))

        name='gradTheta_dwind'
        gradSST_dwind=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))
         
        name='divU10'
        divU10=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))
        
        size=20
        parameters = {'axes.labelsize': size, 'axes.titlesize': size,"axes.labelsize" : size,"xtick.labelsize" : size,\
                "ytick.labelsize" : size,"ytick.major.size" : 5,"ytick.major.width" : 2,"ytick.minor.size" : 3,"ytick.minor.width" : 1,\
       "xtick.major.size" : 5,"xtick.major.width" : 2,"xtick.minor.size" : 3,"xtick.minor.width" : 1}
        
        plt.rcParams.update(parameters)

        fig=plt.figure(figsize=(30,22))

        gs=gridspec.GridSpec(3,2,height_ratios= [2.5,1,1])
        gs.update(bottom=0.03, top=0.905, left=0.08, right=0.92,
                        wspace=0.5, hspace=0.1)
        
        H0=H.isel(depth=slice(0,-1))
        lat0,depth0=H0.lat.data,H0.depth.data
        latt0,depthh0=np.meshgrid(lat0,depth0)
        lat,depth=H.lat.data,H.depth.data
        latt,depthh=np.meshgrid(lat,depth)

        quiverscale=0.6e-4
        quiverwidth=0.0025
        alpha=2

        tx,ty=0.03,0.92
        dcolor=0.013
        labelpad=10
        hscale=1e-3
        
        # --------------------------------------------------------------
        axx=plt.subplot(gs[0,1])
        axx.yaxis.set_label_position('right')
        axx.yaxis.set_ticks_position('right')

        axx.text(0.5,1.18,'2-month mean',weight="bold",size=25,horizontalalignment='center',verticalalignment='center',\
                transform = axx.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
        axx.text(tx,ty+0.05,'c',weight="bold",size=25,horizontalalignment='center',verticalalignment='center',\
                transform = axx.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
        xfront=36.5
        dx=150
        conv=(np.pi*Rt*1e-3)/180
        lat1,lat2=xfront-0.5*dx/conv,xfront+0.5*dx/conv
        ypos=60*hscale
        arr=axx.annotate('',size=15,xy=(lat1, ypos), xytext=(lat2, ypos),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='b', lw=4),annotation_clip=False)
        axx.text((lat2+lat1)/2,ypos+0.01,str(dx)+' km',size=20,color='b',horizontalalignment='center',verticalalignment='bottom')
        
        axx.scatter(xfront,0,marker='*',color='b',s=400,clip_on=False)
        
        
        xfront=39.1
        dx=100
        conv=(np.pi*Rt*1e-3)/180
        lat1,lat2=xfront-0.5*dx/conv,xfront+0.5*dx/conv
        ypos=60*hscale
        arr=axx.annotate('',size=15,xy=(lat1, ypos), xytext=(lat2, ypos),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='b', lw=4),annotation_clip=False)
        axx.text((lat2+lat1)/2,ypos+0.01,str(dx)+' km',size=20,color='b',horizontalalignment='center',verticalalignment='bottom')
        
        frontloc=axx.scatter(xfront,0,marker='*',color='k',s=400,clip_on=False)
        axx.scatter(xfront,0,marker='*',color='b',s=400,clip_on=False)
        
        
        c3=axx.contour(latt0,H0.H.data*hscale,
                        (DTHDT_DELP.DTHDT_DELP*3600*24).data,\
                levels=[12,13,15,20],colors='k',alpha=1,extend='max',linewidths=3)
        axx.clabel(c3, inline=1, fontsize=20)
        leg3,_=c3.legend_elements()#c3.set_label('Dabatic heating')
        
        c3=axx.contourf(latt0,H0.H.data*hscale,
                        ((wrho.wrho-wrho.wrho.mean(dim='lat').data[:,None])*3600*24).data,\
                levels=np.arange(-1000,1001,50),cmap=cmocean.cm.balance,alpha=1,extend='both',linewidths=3) #cmocean.cm.thermal

        c4=axx.contour(latt0,H0.H.data*hscale,
                        (-divV.V*QV.QV*1e3*3600*24).data,\
                levels=[5,7,9],colors='r',alpha=1,extend='max',linewidths=3)
        axx.clabel(c4, inline=1, fontsize=20)
        leg4,_=c4.legend_elements()
        
        tick_font_size=20
        pos1 = axx.get_position()
        cax = fig.add_axes([pos1.x0, pos1.y1+dcolor-0.004, pos1.x1-pos1.x0, +dcolor-0.004])
        cbar=plt.colorbar(mappable=c3,cax=cax,orientation="horizontal",ticklocation='top',aspect=35)
        # cbar.set_label("Pressure weighted dTHdt due to moist [Pa.K.day$^{-1}$]",labelpad=labelpad)
        cbar.set_label("Total vertical mass flux anomaly $(\\rho w)_{anom}$ [kg.m$^{-2}$.day$^{-1}$]",labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)
        
        v=V.V.where(H.H<zmax).data[:,::alpha]/((latmax-latmin)*110e3)
        w=W.W.where(H.H<zmax).data[:,::alpha]/zmax
        vmean=V.V.where(H.H<zmax).mean(dim='lat').data[:,None]/((latmax-latmin)*110e3)
        wmean=W.W.where(H.H<zmax).mean(dim='lat').data[:,None]/zmax
        
        print(np.nanmax(v-vmean)*((latmax-latmin)*110e3),np.nanmax(w-wmean)*zmax)
        q=axx.quiver(latt0[:,::alpha],H0.H.data[:,::alpha]*hscale,\
                v-vmean,
                w-wmean,zorder=2,\
                width=quiverwidth,scale=quiverscale)
        
        lkey=1/((latmax-latmin)*110e3)
        txq,tyq=0.855,0.915
        tt=axx.text(txq+0.038,tyq+0.02,'                \n\n',weight="bold",size=23,horizontalalignment='center',verticalalignment='center',\
                transform = axx.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        tt.set_zorder(2)
        
        txq,tyq=0.855,0.97
        tt=axx.text(txq+0.038,tyq,'Wind anomaly',weight="bold",size=15,horizontalalignment='center',verticalalignment='center',\
                transform = axx.transAxes)
        tt.set_zorder(2)
        qt=axx.quiverkey(q, X=txq, Y=tyq-0.07, U=lkey,
             label=str(lkey*((latmax-latmin)*110e3))+' m.s$^{-1}$', labelpos='E',fontproperties={'weight': 'bold','size': 15},zorder=3)#,bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        t = qt.text.set_zorder(3)
        
        lkey=0.01/zmax
        qt=axx.quiverkey(q, X=txq-0.04, Y=tyq-0.03, U=lkey,
             label=str(lkey*zmax*1e2)+' cm.s$^{-1}$', angle=90,labelsep=0.1,labelpos='E',fontproperties={'weight': 'bold','size': 15},zorder=3)#,bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        t = qt.text.set_zorder(3)
        
        pblh,=axx.plot(lat,PBLH.PBLH.data*hscale,c='purple',linewidth=5,linestyle="--")     
        axx.set_ylim(0,zmax*hscale)
        
        axx.set_ylabel('Height [km]',labelpad=15)

        axx.xaxis.set_ticks([36,37,38,39,40]) 
        labels = [item.get_text()+'°N' for item in axx.get_xticklabels()]
        axx.set_xticklabels(labels) 
        

        
        # --------------------------------------------------------------
        axx1= plt.subplot(gs[1,1],sharex=axx)
        axx1.text(tx,ty,'d',weight="bold",size=25,horizontalalignment='center',verticalalignment='center',\
                transform = axx1.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
        xfront=36.5
        dx=150
        conv=(np.pi*Rt*1e-3)/180
        lat1,lat2=xfront-0.5*dx/conv,xfront+0.5*dx/conv
        ypos=10.8
        
        arr=axx1.annotate('',size=15,xy=(lat1, ypos), xytext=(lat2, ypos),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='b', lw=4),annotation_clip=False)
        axx1.text((lat2+lat1)/2,ypos+0.1,str(dx)+' km',size=20,color='b',horizontalalignment='center',verticalalignment='bottom')
        
        axx1.scatter(xfront,10.5,marker='*',color='b',s=400,clip_on=False)
        
        xfront=39.1
        dx=100
        conv=(np.pi*Rt*1e-3)/180
        lat1,lat2=xfront-0.5*dx/conv,xfront+0.5*dx/conv
        ypos=10.8
        arr=axx1.annotate('',size=15,xy=(lat1, ypos), xytext=(lat2, ypos),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='b', lw=4),annotation_clip=False)
        axx1.text((lat2+lat1)/2,ypos+0.1,str(dx)+' km',color='b',size=20,horizontalalignment='center',verticalalignment='bottom')
        
        axx1.scatter(xfront,10.5,marker='*',color='b',s=400,clip_on=False)
        
        
        axx2 = axx1.twinx()
        axx2.spines.left.set_position(("axes", -0.11))
        axx2.yaxis.set_label_position('left')
        axx2.yaxis.set_ticks_position('left')
        axx2__ = axx1.twinx()
        axx2_ = axx1.twinx()
        axx2_.spines.right.set_position(("axes", 1.13))
        linewidth=5
        SST['Theta'].plot(ax=axx1,c='k',linewidth=linewidth,label='SST')
        EFLUX.EFLUX.plot(ax=axx2,c='r',linewidth=linewidth,label='EFLUX',linestyle='dashed')
        CAPE.CAPE.plot(ax=axx2__,c='g',linewidth=linewidth,label='CAPE')
        (PRECTOT.PRECTOT*24*3600).plot(ax=axx2_,c='b',linewidth=linewidth,label='Convective precipitation')

        axx1.set_title('')
        axx2.set_title('')
        axx1.set_xlabel('')
        axx1.set_ylim(10.5,19.5)
        axx1.set_ylabel('SST [°C]')
        
        axx2.set_ylabel('LHF [W.m$^{-2}$]')
        axx2.set_ylim(100,260)
        axx2.yaxis.label.set_color(axx2.get_lines()[0].get_color())
        axx2.tick_params(axis='y', colors=axx2.get_lines()[0].get_color())        
        
        axx1.yaxis.label.set_color(axx1.get_lines()[0].get_color())
        axx1.tick_params(axis='y', colors=axx1.get_lines()[0].get_color())        
        axx2__.set_title('')
        axx2__.set_ylabel('CAPE [J.kg$^{-1}$]')
        axx2__.yaxis.label.set_color(axx2__.get_lines()[0].get_color())
        axx2__.tick_params(axis='y', colors=axx2__.get_lines()[0].get_color())
        
        axx2_.set_title('')
        axx2_.set_ylabel('Tot. precipitation [mm.day$^{-1}$]')
        axx2_.yaxis.label.set_color(axx2_.get_lines()[0].get_color())
        axx2_.tick_params(axis='y', colors=axx2_.get_lines()[0].get_color())
        
        #----------------------------------------------------------------------------------------------
        axx3= plt.subplot(gs[2,1],sharex=axx)
        axx3.text(tx,ty,'c',weight="bold",size=25,horizontalalignment='center',verticalalignment='center',\
                transform = axx3.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
        xfront=36.5
        dx=150
        conv=(np.pi*Rt*1e-3)/180
        lat1,lat2=xfront-0.5*dx/conv,xfront+0.5*dx/conv
        ypos=-6.8
        arr=axx3.annotate('',size=15,xy=(lat1, ypos), xytext=(lat2, ypos),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='b', lw=4),annotation_clip=False)
        axx3.text((lat2+lat1)/2,ypos+0.01,str(dx)+' km',color='b',size=20,horizontalalignment='center',verticalalignment='bottom')
        
        axx3.scatter(xfront,-7.2,marker='*',color='b',s=400,clip_on=False)
        
        xfront=39.1
        dx=100
        conv=(np.pi*Rt*1e-3)/180
        lat1,lat2=xfront-0.5*dx/conv,xfront+0.5*dx/conv
        ypos=-6.8
        arr=axx3.annotate('',size=15,xy=(lat1, ypos), xytext=(lat2, ypos),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='b', lw=4),annotation_clip=False)
        axx3.text((lat2+lat1)/2,ypos+0.1,str(dx)+' km',color='b',size=20,horizontalalignment='center',verticalalignment='bottom')
        
        axx3.scatter(xfront,-7.2,marker='*',color='b',s=400,clip_on=False)
        
        axx4 = axx3.twinx()
        axx4_= axx3.twinx()
        axx4_.spines.right.set_position(("axes", 1.13))

        p1=(gradSST_dwind.gradTheta_dwind*1e5)\
                .plot(ax=axx3,c='k',linewidth=linewidth,label='$\\nabla$SST downwind')
        p2=(W_2000mean.W*1e2).plot(ax=axx4_,c='b',linestyle='solid',linewidth=linewidth,label='<W(z<2000m)>')
        p4=(divU10.divU10*1e5).plot(ax=axx4,c='r',linestyle='dashed',linewidth=linewidth,label='$\\nabla$U10')

        axx3.set_title('')
        axx3.set_ylabel('$\mathbf{k}.\mathbf{\\nabla}$SST [°C.(100km)$^{-1}$]')
        axx3.set_ylim(-7.2,7.2)
        axx4.set_ylabel('$\mathbf{\\nabla}.\mathbf{u_{10}}$ [m.s$^{-1}$.(100km)$^{-1}$]')
        axx4.set_title('')
        axx4.set_ylim(-4.4,4.4)
        axx3.set_xlabel('')
        axx4_.set_ylim(-2,2)
        axx4_.set_ylabel('$<w(z<$2km$)>$ [cm.s$^{-1}$]')
        axx4_.set_title('')
        axx4_.yaxis.label.set_color(axx4_.get_lines()[0].get_color())
        axx4_.tick_params(axis='y', colors=axx4_.get_lines()[0].get_color())
        axx3.yaxis.label.set_color(axx3.get_lines()[0].get_color())
        axx3.tick_params(axis='y', colors=axx3.get_lines()[0].get_color())
        axx4.yaxis.label.set_color(axx4.get_lines()[0].get_color())
        axx4.tick_params(axis='y', colors=axx4.get_lines()[0].get_color())
        
        
        #################### Cold sector ##################################
        collection=collection2
        
        latmin,latmax=40,41.5
        rollat=2
        rolz=1
        Rt=6356.7523e3
        
        name='H'
        H=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))
        
        name='wrho'
        wrho=xr.open_dataset(folder+name+collection).rolling(lat=rollat,depth=rolz,center=True).mean().sel(lat=slice(latmin,latmax))
        
        name='DTHDT_DELP'
        DTHDT_DELP=xr.open_dataset(folder+name+collection).rolling(lat=rollat,depth=rolz,center=True).mean().sel(lat=slice(latmin,latmax))
        
        name='P'
        P=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))
        
        name='PRECTOT'
        PRECTOT=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))

        name='QV'
        QV=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))

        name='W'
        W=xr.open_dataset(folder+name+collection).rolling(lat=rollat,depth=rolz,center=True).mean().sel(lat=slice(latmin,latmax))

        name='V'
        V=xr.open_dataset(folder+name+collection).rolling(lat=rollat,depth=rolz,center=True).mean().sel(lat=slice(latmin,latmax))

        divV=V.differentiate("lat")*180/(np.pi*Rt)
        
        name='Theta_ocean'
        SST=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))
       
        name='EFLUX'
        EFLUX=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))

        name='CAPE'
        CAPE=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))

        name='W_2000mean'
        W_2000mean=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))

        name='PBLH'
        PBLH=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))

        name='gradTheta_dwind'
        gradSST_dwind=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))
         
        name='divU10'
        divU10=xr.open_dataset(folder+name+collection).sel(lat=slice(latmin,latmax))

        H0=H.isel(depth=slice(0,-1))
        lat0,depth0=H0.lat.data,H0.depth.data
        latt0,depthh0=np.meshgrid(lat0,depth0)
        lat,depth=H.lat.data,H.depth.data
        latt,depthh=np.meshgrid(lat,depth)

        quiverscale=1.5e-4
        quiverwidth=0.0028
        alpha=1
        tx,ty=0.03,0.92
        dcolor=0.013
        labelpad=10
        
        # --------------------------------------------------------------
        axx=plt.subplot(gs[0,0])

        axx.annotate('',size=15,xy=(41.8, 5000*hscale), xytext=(41.8, -5000*hscale),ha='center', va='top',
            arrowprops=dict(arrowstyle='-',color='k', lw=4),annotation_clip=False)
        
        
        axx.text(0.5,1.18,'5-day mean',weight="bold",size=25,horizontalalignment='center',verticalalignment='center',\
                transform = axx.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))

        axx.text(tx,ty+0.05,'a',weight="bold",size=25,horizontalalignment='center',verticalalignment='center',\
                transform = axx.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
        xfront=40.65
        dx=50
        conv=(np.pi*Rt*1e-3)/180
        lat1,lat2=xfront-0.5*dx/conv,xfront+0.5*dx/conv
        ypos=60*hscale
        arr=axx.annotate('',size=15,xy=(lat1, ypos), xytext=(lat2, ypos),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='r', lw=4),annotation_clip=False)
        axx.text((lat2+lat1)/2,ypos+0.01,str(dx)+' km',color='r',size=20,horizontalalignment='center',verticalalignment='bottom')
        
        axx.scatter(xfront,0,marker='*',color='r',s=400,clip_on=False)
        
        c3=axx.contour(latt0,H0.H.data*hscale,
                        (DTHDT_DELP.DTHDT_DELP*3600*24).data,\
                levels=[16,18,20,23],colors='k',alpha=1,extend='max',linewidths=3)#[12,13,15,20]
        axx.clabel(c3, inline=1, fontsize=20)
        
        c3=axx.contourf(latt0,H0.H.data*hscale,
                        ((wrho.wrho-wrho.wrho.mean(dim='lat').data[:,None])*3600*24).data,\
                levels=np.arange(-3200,3201,100),cmap=cmocean.cm.balance,alpha=1,extend='both',linewidths=3)#np.arange(-1000,1001,50)
        
        c4=axx.contour(latt0,H0.H.data*hscale,
                        (-divV.V*QV.QV*1e3*3600*24).data,\
                levels=[3,5],colors='r',alpha=1,extend='max',linewidths=4)
        axx.clabel(c4, inline=1, fontsize=20)
        
        tick_font_size=20
        
        pos1 = axx.get_position()
        cax = fig.add_axes([pos1.x0, pos1.y1+dcolor-0.004, pos1.x1-pos1.x0, +dcolor-0.004])
        cbar=plt.colorbar(mappable=c3,cax=cax,ticks=[-3000,-2000,-1000,0,1000,2000,3000],orientation="horizontal",ticklocation='top',aspect=35)
        cbar.set_label("Total vertical mass flux anomaly $(\\rho w)_{anom}$ [kg.m$^{-2}$.day$^{-1}$]",labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)
        
        v=V.V.where(H.H<zmax).data[:,::alpha]/((latmax-latmin)*110e3)
        w=W.W.where(H.H<zmax).data[:,::alpha]/zmax
        vmean=V.V.where(H.H<zmax).mean(dim='lat').data[:,None]/((latmax-latmin)*110e3)
        wmean=W.W.where(H.H<zmax).mean(dim='lat').data[:,None]/zmax
        print(np.nanmax(v-vmean)*((latmax-latmin)*110e3),np.nanmax(w-wmean)*zmax)
        q=axx.quiver(latt0[:,::alpha],H0.H.data[:,::alpha]*hscale,\
                v-vmean,
                w-wmean,zorder=2,\
                width=quiverwidth,scale=quiverscale)
        
        lkey=1/((latmax-latmin)*110e3)
        txq,tyq=0.855,0.915
        tt=axx.text(txq+0.038,tyq+0.02,'                \n\n',weight="bold",size=23,horizontalalignment='center',verticalalignment='center',\
                transform = axx.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        tt.set_zorder(2)
        
        txq,tyq=0.855,0.97
        tt=axx.text(txq+0.038,tyq,'Wind anomaly',weight="bold",size=15,horizontalalignment='center',verticalalignment='center',\
                transform = axx.transAxes)
        tt.set_zorder(2)
        qt=axx.quiverkey(q, X=txq, Y=tyq-0.07, U=lkey,
             label=str(lkey*((latmax-latmin)*110e3))+' m.s$^{-1}$', labelpos='E',fontproperties={'weight': 'bold','size': 15},zorder=3)#,bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        t = qt.text.set_zorder(3)
        
        lkey=0.02/zmax
        qt=axx.quiverkey(q, X=txq-0.04, Y=tyq-0.03, U=lkey,
             label=str(lkey*zmax*1e2)+' cm.s$^{-1}$', angle=90,labelsep=0.1,labelpos='E',fontproperties={'weight': 'bold','size': 15},zorder=3)#,bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        t = qt.text.set_zorder(3)
        
        axx.plot(lat,PBLH.PBLH.data*hscale,c='purple',linewidth=5,linestyle="--")     
        axx.set_ylim(0,zmax*hscale)
        axx.set_ylabel('Height [km]',labelpad=15)

        axx.xaxis.set_ticks([40.2,40.4,40.6,40.8,41,41.2]) 
        labels = [item.get_text()+'°N' for item in axx.get_xticklabels()]
        axx.set_xticklabels(labels) 
        
        arrow=axx.scatter([], [], marker='$\longleftrightarrow$', c="black", s=1500)
        legend=axx.legend([leg3[0],pblh,leg4[0],arrow,frontloc],\
                ['Tot. diabatic heating','PBLH','Moisture convergence','SST front size','SST front location'],\
                 bbox_to_anchor=(1.25, 0.5), loc='center',fancybox=True,shadow=True,fontsize=22,edgecolor='black')
        legend.get_frame().set_alpha(None)
        
        # --------------------------------------------------------------
        axx1= plt.subplot(gs[1,0],sharex=axx)

        l=axx1.text(tx,ty,'b',weight="bold",size=25,horizontalalignment='center',verticalalignment='center',\
                transform = axx1.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))#.set_zorder(5)
        l.set_zorder(3)
        
        xfront=40.65
        dx=50
        conv=(np.pi*Rt*1e-3)/180
        lat1,lat2=xfront-0.5*dx/conv,xfront+0.5*dx/conv
        ypos=5.3
        
        arr=axx1.annotate('',size=15,xy=(lat1, ypos), xytext=(lat2, ypos),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='r', lw=4),annotation_clip=False)
        axx1.text((lat2+lat1)/2,ypos+0.1,str(dx)+' km',size=20,color='r',horizontalalignment='center',verticalalignment='bottom')
        
        axx1.scatter(xfront,5,marker='*',color='r',s=400,clip_on=False)
        
        axx2 = axx1.twinx()
        axx2.spines.left.set_position(("axes", -0.11))
        axx2.yaxis.set_label_position('left')
        axx2.yaxis.set_ticks_position('left')
        axx2__ = axx1.twinx()
        axx2_ = axx1.twinx()
        axx2_.spines.right.set_position(("axes", 1.13))
        linewidth=5
        SST['Theta'].plot(ax=axx1,c='k',linewidth=linewidth,label='SST',zorder=2)
        EFLUX.EFLUX.plot(ax=axx2,c='r',linewidth=linewidth,label='EFLUX',linestyle='dashed',zorder=1)
        CAPE.CAPE.plot(ax=axx2__,c='g',linewidth=linewidth,label='CAPE',zorder=2)
        (PRECTOT.PRECTOT*24*3600).plot(ax=axx2_,c='b',linewidth=linewidth,label='PRECTOT',zorder=2)

        axx1.set_title('')
        axx2.set_title('')
        axx1.set_xlabel('')
        axx1.set_ylabel('SST [°C]')
        axx1.set_ylim(5,14)
        axx2.set_ylabel('LHF [W.m$^{-2}$]')
        axx2.set_ylim(0,185)
        axx2.yaxis.label.set_color(axx2.get_lines()[0].get_color())
        axx2.tick_params(axis='y', colors=axx2.get_lines()[0].get_color())        
        
        axx1.yaxis.label.set_color(axx1.get_lines()[0].get_color())
        axx1.tick_params(axis='y', colors=axx1.get_lines()[0].get_color())        
        axx2__.set_title('')
        axx2__.set_ylabel('CAPE [J.kg$^{-1}$]')
        axx2__.yaxis.label.set_color(axx2__.get_lines()[0].get_color())
        axx2__.tick_params(axis='y', colors=axx2__.get_lines()[0].get_color())
        
        axx2_.set_title('')
        axx2_.set_ylabel('Tot. precipitation [mm.day$^{-1}$]')
        axx2_.yaxis.label.set_color(axx2_.get_lines()[0].get_color())
        axx2_.tick_params(axis='y', colors=axx2_.get_lines()[0].get_color())
        
        # --------------------------------------------------------------
        axx3= plt.subplot(gs[2,0],sharex=axx)
        axx3.text(tx,ty,'f',weight="bold",size=25,horizontalalignment='center',verticalalignment='center',\
                transform = axx3.transAxes,\
                bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        
        xfront=40.65
        dx=50
        conv=(np.pi*Rt*1e-3)/180
        lat1,lat2=xfront-0.5*dx/conv,xfront+0.5*dx/conv
        ypos=-15
        
        arr=axx3.annotate('',size=15,xy=(lat1, ypos), xytext=(lat2, ypos),ha='center', va='center',
            arrowprops=dict(arrowstyle='<->',color='r', lw=4),annotation_clip=False)
        axx3.text((lat2+lat1)/2,ypos+0.1,str(dx)+' km',color='r',size=20,horizontalalignment='center',verticalalignment='bottom')
        
        axx3.scatter(xfront,-16,marker='*',color='r',s=400,clip_on=False)
        
        axx4 = axx3.twinx()
        axx4_= axx3.twinx()
        axx4_.spines.right.set_position(("axes", 1.13))

        p1=(gradSST_dwind.gradTheta_dwind*1e5)\
                .plot(ax=axx3,c='k',linewidth=linewidth,label='$\\nabla$SST downwind')
        p2=(W_2000mean.W*1e2).plot(ax=axx4_,c='b',linestyle='solid',linewidth=linewidth,label='<W(z<2000m)>')
        p4=(divU10.divU10*1e5).plot(ax=axx4,c='r',linestyle='dashed',linewidth=linewidth,label='$\\nabla$U10')

        
        axx3.set_title('')
        axx3.set_ylabel('$\mathbf{k}.\mathbf{\\nabla}$SST [°C.(100km)$^{-1}$]')
        axx3.set_ylim(-16,16)
        axx4.set_ylabel('$\mathbf{\\nabla}.\mathbf{u_{10}}$ [m.s$^{-1}$.(100km)$^{-1}$]')
        axx4.set_title('')
        axx4.set_ylim(-8.5,8.5)
        axx4.set_yticks([-8,-4,0,4,8])
        axx3.set_xlabel('')
        axx4_.set_ylim(-2.3,2.3)
        axx4_.set_ylabel('$<w(z<$2km$)>$ [cm.s$^{-1}$]')
        axx4_.set_title('')
        axx4_.yaxis.label.set_color(axx4_.get_lines()[0].get_color())
        axx4_.tick_params(axis='y', colors=axx4_.get_lines()[0].get_color())
        axx3.yaxis.label.set_color(axx3.get_lines()[0].get_color())
        axx3.tick_params(axis='y', colors=axx3.get_lines()[0].get_color())
        axx4.yaxis.label.set_color(axx4.get_lines()[0].get_color())
        axx4.tick_params(axis='y', colors=axx4.get_lines()[0].get_color())

       
lonsel=152.5
folder='/path/to/your/extracted_data/2Dv/'
collection='_2Dv_decjan_NWsup0_fullsection.nc'
collection2='_2Dv_16dec21dec_NWsup0_fullsection.nc'
plott('title',folder,lonsel,collection,collection2)
plt.savefig('/path/to/your/plots/ExtDataFig6.png')
plt.close()
        

