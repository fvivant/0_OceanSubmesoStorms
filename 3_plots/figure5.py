# Coded by Felix Vivant (Scripps-UCSD and ENS Paris-Saclay)
# Program to plot Fig.5

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec


# --------------- Plot function
def plott(folder, collection, collection2, figsize, scaling):

        zmax = 4000
        latmin, latmax = 40, 41.5
        rollat = 2
        rolz = 1
        Rt = 6356.7523e3

        # -------- Open data
        name = 'H'
        H = xr.open_dataset(folder + name +
                            collection).sel(lat=slice(latmin, latmax))

        name = 'DTHDTCN'
        DTHDTCN = xr.open_dataset(folder + name +
                                  collection).sel(lat=slice(latmin, latmax))

        name = 'CNVMFC'
        CNVMFC = xr.open_dataset(folder + name +
                                 collection).sel(lat=slice(latmin, latmax))

        name = 'QV'
        QV = xr.open_dataset(folder + name +
                             collection).sel(lat=slice(latmin, latmax))

        name = 'W'
        W = xr.open_dataset(folder + name + collection).rolling(
            lat=rollat, depth=rolz,
            center=True).mean().sel(lat=slice(latmin, latmax))

        name = 'V'
        V = xr.open_dataset(folder + name + collection).rolling(
            lat=rollat, depth=rolz,
            center=True).mean().sel(lat=slice(latmin, latmax))

        name = 'Theta_ocean'
        SST = xr.open_dataset(folder + name +
                              collection).sel(lat=slice(latmin, latmax))

        name = 'EFLUX'
        EFLUX = xr.open_dataset(folder + name +
                                collection).sel(lat=slice(latmin, latmax))

        name = 'PRECCON'
        PRECCON = xr.open_dataset(folder + name +
                                  collection).sel(lat=slice(latmin, latmax))

        name = 'CAPE'
        CAPE = xr.open_dataset(folder + name +
                               collection).sel(lat=slice(latmin, latmax))

        name = 'W_2000mean'
        W_2000mean = xr.open_dataset(folder + name +
                                     collection).sel(lat=slice(latmin, latmax))

        name = 'PBLH'
        PBLH = xr.open_dataset(folder + name +
                               collection).sel(lat=slice(latmin, latmax))

        name = 'gradTheta_dwind'
        gradSST_dwind = xr.open_dataset(folder + name + collection).sel(
            lat=slice(latmin, latmax))

        name = 'divU10'
        divU10 = xr.open_dataset(folder + name +
                                 collection).sel(lat=slice(latmin, latmax))

        # -------- Defintion of the plot
        size = 6 * scaling
        largesize = 8 * scaling
        smallsize = 5 * scaling
        parameters = {
            'axes.labelsize': size,
            'axes.titlesize': size,
            "axes.labelsize": size,
            "xtick.labelsize": size,
            "ytick.labelsize": size,
            "ytick.major.size": 5,
            "ytick.major.width": 2,
            "ytick.minor.size": 3,
            "ytick.minor.width": 1,
            "xtick.major.size": 5,
            "xtick.major.width": 2,
            "xtick.minor.size": 3,
            "xtick.minor.width": 1
        }

        plt.rcParams.update(parameters)

        fig = plt.figure(figsize=figsize)

        gs = gridspec.GridSpec(3, 2, height_ratios=[2.5, 1, 1])
        gs.update(bottom=0.03,
                  top=0.905,
                  left=0.08,
                  right=0.92,
                  wspace=0.52,
                  hspace=0.12)

        H0 = H.isel(depth=slice(0, -1))
        lat0, depth0 = H0.lat.data, H0.depth.data
        latt0, depthh0 = np.meshgrid(lat0, depth0)
        lat, depth = H.lat.data, H.depth.data
        latt, depthh = np.meshgrid(lat, depth)

        # Parameters
        quiverscale = 1.5e-4
        quiverwidth = 0.0028
        alpha = 1
        tx, ty = 0.03, 0.92
        dcolor = 0.013
        labelpad = 12
        hscale = 1e-3

        # --------------------------------------------------------------
        axx = plt.subplot(gs[0, 0])
        # axx.yaxis.set_label_position('right')
        # axx.yaxis.set_ticks_position('right')

        axx.text(0.5,
                 1.2,
                 'Warm sector',
                 weight="bold",
                 size=largesize,
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform=axx.transAxes,
                 bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))

        axx.annotate('',
                     size=15,
                     xy=(41.8, 5000 * hscale),
                     xytext=(41.8, -4100 * hscale),
                     ha='center',
                     va='top',
                     arrowprops=dict(arrowstyle='-', color='k', lw=4),
                     annotation_clip=False)

        axx.text(tx + 0.01,
                 ty + 0.03,
                 'a',
                 weight="bold",
                 size=largesize,
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform=axx.transAxes,
                 bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))

        xfront = 40.65
        dx = 50
        conv = (np.pi * Rt * 1e-3) / 180
        lat1, lat2 = xfront - 0.5 * dx / conv, xfront + 0.5 * dx / conv
        ypos = 60 * hscale
        col = 'r'
        arr = axx.annotate('',
                           size=15,
                           xy=(lat1, ypos),
                           xytext=(lat2, ypos),
                           ha='center',
                           va='center',
                           arrowprops=dict(arrowstyle='<->', color=col, lw=4),
                           annotation_clip=False)
        axx.text((lat2 + lat1) / 2,
                 ypos + 0.01,
                 str(dx) + ' km',
                 size=smallsize,
                 color=col,
                 horizontalalignment='center',
                 verticalalignment='bottom')

        frontloc = axx.scatter(xfront,
                               0,
                               marker='*',
                               color=col,
                               s=400,
                               clip_on=False)
        axx.scatter(xfront, 0, marker='*', color=col, s=400, clip_on=False)

        c3=axx.contour(latt0,H0.H.data*hscale,
                        (DTHDTCN.DTHDTCN*3600*24).data,\
                levels=[8,12,15],colors='k',alpha=1,extend='max',linewidths=4)
        axx.clabel(c3, inline=1, fontsize=20)
        leg3, _ = c3.legend_elements()

        c3=axx.contourf(latt0,H0.H.data*hscale,
                        (CNVMFC.CNVMFC*3600*24).data,\
                levels=np.arange(0,3501,300),cmap='YlOrRd',alpha=0.9,extend='max') #.where((CNVMFC.CNVMFC*3600*24)>=300)

        tick_font_size = size
        pos1 = axx.get_position()
        cax = fig.add_axes([
            pos1.x0, pos1.y1 + dcolor - 0.004, pos1.x1 - pos1.x0,
            +dcolor - 0.004
        ])
        cbar = plt.colorbar(mappable=c3,
                            cax=cax,
                            orientation="horizontal",
                            ticklocation='top',
                            aspect=35)
        cbar.set_label("Convective Mass Flux [kg.m$^{-2}$.day$^{-1}$]",
                       labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)

        v = V.V.where(H.H < zmax).data[:, ::alpha] / ((latmax - latmin) * 110e3)
        w = W.W.where(H.H < zmax).data[:, ::alpha] / zmax
        vmean = V.V.where(H.H < zmax).mean(dim='lat').data[:, None] / (
            (latmax - latmin) * 110e3)
        wmean = W.W.where(H.H < zmax).mean(dim='lat').data[:, None] / zmax

        q=axx.quiver(latt0[:,::alpha],H0.H.data[:,::alpha]*hscale,\
                v-vmean,
                w-wmean,zorder=2,\
                width=quiverwidth,scale=quiverscale)

        lkey = 1 / ((latmax - latmin) * 110e3)
        txq, tyq = 0.83, 0.895
        tt = axx.text(txq + 0.038,
                      tyq + 0.02,
                      '                 \n\n',
                      weight="bold",
                      size=29,
                      horizontalalignment='center',
                      verticalalignment='center',
                      transform=axx.transAxes,
                      bbox=dict(boxstyle="round",
                                fc="white",
                                ec="black",
                                pad=0.2))
        tt.set_zorder(2)

        # txq, tyq = 0.855, 0.97
        tt = axx.text(txq + 0.04,
                      tyq + 0.065,
                      'Wind anomaly',
                      weight="bold",
                      size=smallsize,
                      horizontalalignment='center',
                      verticalalignment='center',
                      transform=axx.transAxes)
        tt.set_zorder(2)
        qt = axx.quiverkey(
            q,
            X=txq - 0.01,
            Y=tyq - 0.015,
            U=lkey,
            label=str(lkey * ((latmax - latmin) * 110e3)) + ' m.s$^{-1}$',
            labelpos='E',
            fontproperties={
                'weight': 'bold',
                'size': smallsize
            },
            zorder=3
        )  #,bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        t = qt.text.set_zorder(3)

        lkey = 0.02 / zmax
        qt = axx.quiverkey(
            q,
            X=txq - 0.05,
            Y=tyq + 0.025,
            U=lkey,
            label=str(lkey * zmax * 1e2) + ' cm.s$^{-1}$',
            angle=90,
            labelsep=0.1,
            labelpos='E',
            fontproperties={
                'weight': 'bold',
                'size': smallsize
            },
            zorder=3
        )  #,bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        t = qt.text.set_zorder(3)

        pblh, = axx.plot(lat,
                         PBLH.PBLH.data * hscale,
                         c='purple',
                         linewidth=5,
                         linestyle="--")
        axx.set_ylim(0, zmax * hscale)
        axx.set_ylabel('Height [km]', labelpad=15)

        axx.xaxis.set_ticks([40.2, 40.4, 40.6, 40.8, 41, 41.2])
        labels = [item.get_text() + '°N' for item in axx.get_xticklabels()]
        axx.set_xticklabels(labels)

        arrow = axx.scatter([], [],
                            marker='$\longleftrightarrow$',
                            c="red",
                            s=1500)
        legend = axx.legend([leg3[0], pblh, arrow, frontloc], [
            'Conv. diabatic heating', 'PBLH', 'SST front size',
            'SST front location'
        ],
                            bbox_to_anchor=(1.26, 0.5),
                            loc='center',
                            fancybox=True,
                            shadow=True,
                            fontsize=smallsize + 2,
                            edgecolor='black')
        legend.get_frame().set_alpha(None)

        # --------------------------------------------------------------
        axx1 = plt.subplot(gs[1, 0], sharex=axx)

        xfront = 40.65
        dx = 50
        conv = (np.pi * Rt * 1e-3) / 180
        lat1, lat2 = xfront - 0.5 * dx / conv, xfront + 0.5 * dx / conv
        ypos = 5.3

        arr = axx1.annotate('',
                            size=15,
                            xy=(lat1, ypos),
                            xytext=(lat2, ypos),
                            ha='center',
                            va='center',
                            arrowprops=dict(arrowstyle='<->', color=col, lw=4),
                            annotation_clip=False)
        axx1.text((lat2 + lat1) / 2,
                  ypos + 0.1,
                  str(dx) + ' km',
                  size=smallsize,
                  color=col,
                  horizontalalignment='center',
                  verticalalignment='bottom')

        axx1.scatter(xfront, 5, marker='*', color=col, s=400, clip_on=False)

        axx2 = axx1.twinx()
        axx2.spines.left.set_position(("axes", -0.11))
        axx2.yaxis.set_label_position('left')
        axx2.yaxis.set_ticks_position('left')
        axx2__ = axx1.twinx()
        axx2_ = axx1.twinx()
        axx2_.spines.right.set_position(("axes", 1.13))
        linewidth = 5
        SST['Theta'].plot(ax=axx1,
                          c='k',
                          linewidth=linewidth,
                          label='SST',
                          zorder=2)
        EFLUX.EFLUX.plot(ax=axx2,
                         c='r',
                         linewidth=linewidth,
                         label='EFLUX',
                         linestyle='dashed',
                         zorder=1)
        CAPE.CAPE.plot(ax=axx2__,
                       c='g',
                       linewidth=linewidth,
                       label='CAPE',
                       zorder=2)
        (PRECCON.PRECCON * 24 * 3600).plot(ax=axx2_,
                                           c='b',
                                           linewidth=linewidth,
                                           label='PRECCON',
                                           zorder=2)

        l = axx2_.text(tx + 0.01,
                       ty - 0.03,
                       'b',
                       weight="bold",
                       size=largesize,
                       horizontalalignment='center',
                       verticalalignment='center',
                       transform=axx1.transAxes,
                       bbox=dict(boxstyle="round",
                                 fc="white",
                                 ec="black",
                                 pad=0.2))  #.set_zorder(5)

        axx1.set_title('')
        axx2.set_title('')
        axx1.set_xlabel('')
        axx1.set_ylabel('SST [°C]')
        axx1.set_ylim(5, 14)
        axx2.set_ylabel('LHF [W.m$^{-2}$]')
        axx2.set_ylim(0, 185)
        axx2.yaxis.label.set_color(axx2.get_lines()[0].get_color())
        axx2.tick_params(axis='y', colors=axx2.get_lines()[0].get_color())

        axx1.yaxis.label.set_color(axx1.get_lines()[0].get_color())
        axx1.tick_params(axis='y', colors=axx1.get_lines()[0].get_color())
        axx2__.set_title('')
        axx2__.set_ylabel('CAPE [J.kg$^{-1}$]')
        axx2__.yaxis.label.set_color(axx2__.get_lines()[0].get_color())
        axx2__.tick_params(axis='y', colors=axx2__.get_lines()[0].get_color())

        axx2_.set_title('')
        axx2_.set_ylim(0, 15)
        axx2_.set_yticks([0, 2, 4, 6, 8, 10, 12, 14])
        axx2_.set_ylabel('Conv. precip. [mm.day$^{-1}$]', fontsize=size)
        axx2_.yaxis.label.set_color(axx2_.get_lines()[0].get_color())
        axx2_.tick_params(axis='y', colors=axx2_.get_lines()[0].get_color())

        # --------------------------------------------------------------
        axx3 = plt.subplot(gs[2, 0], sharex=axx)
        axx3.text(tx + 0.01,
                  ty - 0.03,
                  'c',
                  weight="bold",
                  size=largesize,
                  horizontalalignment='center',
                  verticalalignment='center',
                  transform=axx3.transAxes,
                  bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))

        xfront = 40.65
        dx = 50
        conv = (np.pi * Rt * 1e-3) / 180
        lat1, lat2 = xfront - 0.5 * dx / conv, xfront + 0.5 * dx / conv
        ypos = -15

        arr = axx3.annotate('',
                            size=15,
                            xy=(lat1, ypos),
                            xytext=(lat2, ypos),
                            ha='center',
                            va='center',
                            arrowprops=dict(arrowstyle='<->', color=col, lw=4),
                            annotation_clip=False)
        axx3.text((lat2 + lat1) / 2,
                  ypos + 0.1,
                  str(dx) + ' km',
                  size=smallsize,
                  color=col,
                  horizontalalignment='center',
                  verticalalignment='bottom')

        axx3.scatter(xfront, -16, marker='*', color=col, s=400, clip_on=False)

        axx4 = axx3.twinx()
        axx4_ = axx3.twinx()
        axx4_.spines.right.set_position(("axes", 1.13))

        p1=(gradSST_dwind.gradTheta_dwind*1e5)\
                .plot(ax=axx3,c='k',linewidth=linewidth,label='$\\nabla$SST downwind')
        p2 = (W_2000mean.W * 1e2).plot(ax=axx4_,
                                       c='b',
                                       linestyle='solid',
                                       linewidth=linewidth,
                                       label='<W(z<2000m)>')
        p4 = (divU10.divU10 * 1e5).plot(ax=axx4,
                                        c='r',
                                        linestyle='dashed',
                                        linewidth=linewidth,
                                        label='$\\nabla$U10')

        axx3.set_title('')
        axx3.set_ylabel('$\mathbf{k}.\mathbf{\\nabla}$SST [°C.(100km)$^{-1}$]')
        axx3.set_ylim(-16, 16)
        axx4.set_ylabel(
            '$\mathbf{\\nabla}.\mathbf{u_{10}}$ [m.s$^{-1}$.(100km)$^{-1}$]')
        axx4.set_title('')
        axx4.set_ylim(-8.5, 8.5)
        axx4.set_yticks([-8, -4, 0, 4, 8])
        axx3.set_xlabel('')
        axx4_.set_ylim(-2.3, 2.3)
        axx4_.set_ylabel('$w(z<$2km$)$ [cm.s$^{-1}$]')
        axx4_.set_title('')
        axx4_.yaxis.label.set_color(axx4_.get_lines()[0].get_color())
        axx4_.tick_params(axis='y', colors=axx4_.get_lines()[0].get_color())
        axx3.yaxis.label.set_color(axx3.get_lines()[0].get_color())
        axx3.tick_params(axis='y', colors=axx3.get_lines()[0].get_color())
        axx4.yaxis.label.set_color(axx4.get_lines()[0].get_color())
        axx4.tick_params(axis='y', colors=axx4.get_lines()[0].get_color())

        ################################################################################
        collection = collection2

        zmax = 4000
        latmin, latmax = 40, 41.5
        rollat = 2
        rolz = 1
        Rt = 6356.7523e3

        name = 'H'
        H = xr.open_dataset(folder + name +
                            collection).sel(lat=slice(latmin, latmax))

        name = 'DTHDTCN'
        DTHDTCN = xr.open_dataset(folder + name +
                                  collection).sel(lat=slice(latmin, latmax))

        name = 'CNVMFC'
        CNVMFC = xr.open_dataset(folder + name +
                                 collection).sel(lat=slice(latmin, latmax))

        name = 'W'
        W = xr.open_dataset(folder + name + collection).rolling(
            lat=rollat, depth=rolz,
            center=True).mean().sel(lat=slice(latmin, latmax))

        name = 'V'
        V = xr.open_dataset(folder + name + collection).rolling(
            lat=rollat, depth=rolz,
            center=True).mean().sel(lat=slice(latmin, latmax))

        name = 'Theta_ocean'
        SST = xr.open_dataset(folder + name +
                              collection).sel(lat=slice(latmin, latmax))

        name = 'EFLUX'
        EFLUX = xr.open_dataset(folder + name +
                                collection).sel(lat=slice(latmin, latmax))

        name = 'PRECCON'
        PRECCON = xr.open_dataset(folder + name +
                                  collection).sel(lat=slice(latmin, latmax))

        name = 'CAPE'
        CAPE = xr.open_dataset(folder + name +
                               collection).sel(lat=slice(latmin, latmax))

        name = 'W_2000mean'
        W_2000mean = xr.open_dataset(folder + name +
                                     collection).sel(lat=slice(latmin, latmax))

        name = 'PBLH'
        PBLH = xr.open_dataset(folder + name +
                               collection).sel(lat=slice(latmin, latmax))

        name = 'gradTheta_dwind'
        gradSST_dwind = xr.open_dataset(folder + name + collection).sel(
            lat=slice(latmin, latmax))

        name = 'divU10'
        divU10 = xr.open_dataset(folder + name +
                                 collection).sel(lat=slice(latmin, latmax))

        H0 = H.isel(depth=slice(0, -1))
        lat0, depth0 = H0.lat.data, H0.depth.data
        latt0, depthh0 = np.meshgrid(lat0, depth0)
        lat, depth = H.lat.data, H.depth.data
        latt, depthh = np.meshgrid(lat, depth)

        quiverscale = 1.3e-4
        quiverwidth = 0.003
        alpha = 1

        tx, ty = 0.03, 0.92
        dcolor = 0.013
        labelpad = 10
        # --------------------------------------------------------------
        axx = plt.subplot(gs[0, 1])
        axx.yaxis.set_label_position('right')
        axx.yaxis.set_ticks_position('right')

        axx.text(0.5,
                 1.2,
                 'Cold sector',
                 weight="bold",
                 size=largesize,
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform=axx.transAxes,
                 bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))

        axx.text(tx + 0.01,
                 ty + 0.03,
                 'd',
                 weight="bold",
                 size=largesize,
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform=axx.transAxes,
                 bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))

        xfront = 40.65
        dx = 50
        conv = (np.pi * Rt * 1e-3) / 180
        lat1, lat2 = xfront - 0.5 * dx / conv, xfront + 0.5 * dx / conv
        ypos = 60 * hscale
        col = 'r'
        arr = axx.annotate('',
                           size=15,
                           xy=(lat1, ypos),
                           xytext=(lat2, ypos),
                           ha='center',
                           va='center',
                           arrowprops=dict(arrowstyle='<->', color=col, lw=4),
                           annotation_clip=False)
        axx.text((lat2 + lat1) / 2,
                 ypos + 0.01,
                 str(dx) + ' km',
                 size=smallsize,
                 color=col,
                 horizontalalignment='center',
                 verticalalignment='bottom')

        axx.scatter(xfront, 0, marker='*', color=col, s=400, clip_on=False)

        c3=axx.contour(latt0,H0.H.data*hscale,
                        (DTHDTCN.DTHDTCN*3600*24).data,\
                levels=[1,2],colors='k',alpha=1,extend='max',linewidths=4)
        axx.clabel(c3, inline=1, fontsize=20)


        c3=axx.contourf(latt0,H0.H.data*hscale,
                        (CNVMFC.CNVMFC*3600*24).data,\
                levels=np.arange(0,801,50),cmap='YlOrRd',alpha=0.9,extend='max')

        tick_font_size = size
        pos1 = axx.get_position()
        cax = fig.add_axes([
            pos1.x0, pos1.y1 + dcolor - 0.004, pos1.x1 - pos1.x0,
            +dcolor - 0.004
        ])
        cbar = plt.colorbar(mappable=c3,
                            cax=cax,
                            orientation="horizontal",
                            ticklocation='top',
                            aspect=35)
        cbar.set_label("Convective Mass Flux [kg.m$^{-2}$.day$^{-1}$]",
                       labelpad=labelpad)
        cbar.ax.tick_params(labelsize=tick_font_size)

        v = V.V.where(H.H < zmax).data[:, ::alpha] / ((latmax - latmin) * 110e3)
        w = W.W.where(H.H < zmax).data[:, ::alpha] / zmax
        vmean = V.V.where(H.H < zmax).mean(dim='lat').data[:, None] / (
            (latmax - latmin) * 110e3)
        wmean = W.W.where(H.H < zmax).mean(dim='lat').data[:, None] / zmax

        q=axx.quiver(latt0[:,::alpha],H0.H.data[:,::alpha]*hscale,\
                v-vmean,
                w-wmean,zorder=2,\
                width=quiverwidth,scale=quiverscale)

        lkey = 1 / ((latmax - latmin) * 110e3)

        txq, tyq = 0.83, 0.895
        tt = axx.text(txq + 0.038,
                      tyq + 0.02,
                      '                 \n\n',
                      weight="bold",
                      size=29,
                      horizontalalignment='center',
                      verticalalignment='center',
                      transform=axx.transAxes,
                      bbox=dict(boxstyle="round",
                                fc="white",
                                ec="black",
                                pad=0.2))
        tt.set_zorder(2)

        # txq, tyq = 0.855, 0.97
        tt = axx.text(txq + 0.04,
                      tyq + 0.065,
                      'Wind anomaly',
                      weight="bold",
                      size=smallsize,
                      horizontalalignment='center',
                      verticalalignment='center',
                      transform=axx.transAxes)
        tt.set_zorder(2)
        qt = axx.quiverkey(
            q,
            X=txq - 0.01,
            Y=tyq - 0.015,
            U=lkey,
            label=str(lkey * ((latmax - latmin) * 110e3)) + ' m.s$^{-1}$',
            labelpos='E',
            fontproperties={
                'weight': 'bold',
                'size': smallsize
            },
            zorder=3
        )  #,bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        t = qt.text.set_zorder(3)

        lkey = 0.02 / zmax
        qt = axx.quiverkey(
            q,
            X=txq - 0.05,
            Y=tyq + 0.025,
            U=lkey,
            label=str(lkey * zmax * 1e2) + ' cm.s$^{-1}$',
            angle=90,
            labelsep=0.1,
            labelpos='E',
            fontproperties={
                'weight': 'bold',
                'size': smallsize
            },
            zorder=3
        )  #,bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))
        t = qt.text.set_zorder(3)

        axx.plot(lat,
                 PBLH.PBLH.data * hscale,
                 c='purple',
                 linewidth=5,
                 linestyle="--")
        axx.set_ylim(0, zmax * hscale)
        axx.set_ylabel('Height [km]', labelpad=15)

        axx.xaxis.set_ticks([40.2, 40.4, 40.6, 40.8, 41, 41.2])
        labels = [item.get_text() + '°N' for item in axx.get_xticklabels()]
        axx.set_xticklabels(labels)

        # --------------------------------------------------------------
        axx1 = plt.subplot(gs[1, 1], sharex=axx)

        xfront = 40.65
        dx = 50
        conv = (np.pi * Rt * 1e-3) / 180
        lat1, lat2 = xfront - 0.5 * dx / conv, xfront + 0.5 * dx / conv
        ypos = 5.3

        arr = axx1.annotate('',
                            size=15,
                            xy=(lat1, ypos),
                            xytext=(lat2, ypos),
                            ha='center',
                            va='center',
                            arrowprops=dict(arrowstyle='<->', color=col, lw=4),
                            annotation_clip=False)
        axx1.text((lat2 + lat1) / 2,
                  ypos + 0.1,
                  str(dx) + ' km',
                  size=smallsize,
                  color=col,
                  horizontalalignment='center',
                  verticalalignment='bottom')

        axx1.scatter(xfront, 5, marker='*', color=col, s=400, clip_on=False)

        axx2 = axx1.twinx()
        axx2.spines.left.set_position(("axes", -0.11))
        axx2.yaxis.set_label_position('left')
        axx2.yaxis.set_ticks_position('left')
        axx2__ = axx1.twinx()
        axx2_ = axx1.twinx()
        axx2_.spines.right.set_position(("axes", 1.13))
        linewidth = 5
        SST['Theta'].plot(ax=axx1,
                          c='k',
                          linewidth=linewidth,
                          label='SST',
                          zorder=2)
        EFLUX.EFLUX.plot(ax=axx2,
                         c='r',
                         linewidth=linewidth,
                         label='EFLUX',
                         linestyle='dashed',
                         zorder=1)
        CAPE.CAPE.plot(ax=axx2__,
                       c='g',
                       linewidth=linewidth,
                       label='CAPE',
                       zorder=2)
        (PRECCON.PRECCON * 24 * 3600).plot(ax=axx2_,
                                           c='b',
                                           linewidth=linewidth,
                                           label='PRECCON',
                                           zorder=2)

        l = axx2_.text(tx + 0.01,
                       ty - 0.03,
                       'e',
                       weight="bold",
                       size=largesize,
                       horizontalalignment='center',
                       verticalalignment='center',
                       transform=axx1.transAxes,
                       bbox=dict(boxstyle="round",
                                 fc="white",
                                 ec="black",
                                 pad=0.2))  #.set_zorder(5)

        axx1.set_title('')
        axx2.set_title('')
        axx1.set_xlabel('')
        axx1.set_ylabel('SST [°C]')
        axx1.set_ylim(5, 14)
        axx2.set_ylabel('LHF [W.m$^{-2}$]')
        axx2.set_ylim(100, 270)
        axx2.yaxis.label.set_color(axx2.get_lines()[0].get_color())
        axx2.tick_params(axis='y', colors=axx2.get_lines()[0].get_color())
        axx1.yaxis.label.set_color(axx1.get_lines()[0].get_color())
        axx1.tick_params(axis='y', colors=axx1.get_lines()[0].get_color())
        axx2__.set_title('')
        axx2__.set_ylim(30, 75)
        axx2__.set_ylabel('CAPE [J.kg$^{-1}$]')
        axx2__.yaxis.label.set_color(axx2__.get_lines()[0].get_color())
        axx2__.tick_params(axis='y', colors=axx2__.get_lines()[0].get_color())
        axx2_.set_title('')
        axx2_.set_ylim(-0.15, 1)
        axx2_.set_ylabel('Conv. precip. [mm.day$^{-1}$]')
        axx2_.yaxis.label.set_color(axx2_.get_lines()[0].get_color())
        axx2_.tick_params(axis='y', colors=axx2_.get_lines()[0].get_color())

        # --------------------------------------------------------------
        axx3 = plt.subplot(gs[2, 1], sharex=axx)
        axx3.text(tx + 0.01,
                  ty - 0.03,
                  'f',
                  weight="bold",
                  size=largesize,
                  horizontalalignment='center',
                  verticalalignment='center',
                  transform=axx3.transAxes,
                  bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2))

        xfront = 40.65
        dx = 50
        conv = (np.pi * Rt * 1e-3) / 180
        lat1, lat2 = xfront - 0.5 * dx / conv, xfront + 0.5 * dx / conv
        ypos = -27.5

        arr = axx3.annotate('',
                            size=15,
                            xy=(lat1, ypos),
                            xytext=(lat2, ypos),
                            ha='center',
                            va='center',
                            arrowprops=dict(arrowstyle='<->', color=col, lw=4),
                            annotation_clip=False)
        axx3.text((lat2 + lat1) / 2,
                  ypos + 0.1,
                  str(dx) + ' km',
                  size=smallsize,
                  color=col,
                  horizontalalignment='center',
                  verticalalignment='bottom')

        axx3.scatter(xfront, -30, marker='*', color=col, s=400, clip_on=False)

        axx4 = axx3.twinx()
        axx4_ = axx3.twinx()
        axx4_.spines.right.set_position(("axes", 1.13))

        p1=(gradSST_dwind.gradTheta_dwind*1e5)\
                .plot(ax=axx3,c='k',linewidth=linewidth,label='$\\nabla$SST downwind')
        p2 = (W_2000mean.W * 1e2).plot(ax=axx4_,
                                       c='b',
                                       linestyle='solid',
                                       linewidth=linewidth,
                                       label='<W(z<2000m)>')
        p4 = (divU10.divU10 * 1e5).plot(ax=axx4,
                                        c='r',
                                        linestyle='dashed',
                                        linewidth=linewidth,
                                        label='$\\nabla$U10')

        axx3.set_title('')
        axx3.set_ylabel('$\mathbf{k}.\mathbf{\\nabla}$SST [°C.(100km)$^{-1}$]')
        axx3.set_ylim(-30, 30)
        axx4.set_ylabel(
            '$\mathbf{\\nabla}.\mathbf{u_{10}}$ [m.s$^{-1}$.(100km)$^{-1}$]')
        axx4.set_title('')
        axx4.set_ylim(-4, 4)
        axx3.set_xlabel('')
        axx4_.set_ylim(-2.7, 2.7)
        axx4_.set_ylabel('$w(z<$2km$)$ [cm.s$^{-1}$]')
        axx4_.set_title('')
        axx4_.yaxis.label.set_color(axx4_.get_lines()[0].get_color())
        axx4_.tick_params(axis='y', colors=axx4_.get_lines()[0].get_color())
        axx3.yaxis.label.set_color(axx3.get_lines()[0].get_color())
        axx3.tick_params(axis='y', colors=axx3.get_lines()[0].get_color())
        axx4.yaxis.label.set_color(axx4.get_lines()[0].get_color())
        axx4.tick_params(axis='y', colors=axx4.get_lines()[0].get_color())


# ----------- Plot
lonsel = 152.5
# folder = '/path/to/your/extracted_data/2Dv/'
folder = '/nobackup/fvivant/DATA/KExt/winter_decjanfev/2D_152_5/'
collection = '_2Dv_16dec21dec_NWsup0_fullsection.nc'
collection2 = '_2Dv_16dec21dec_SW_fullsection.nc'

figsize = (30, 22)
fig_width = figsize[0]  # inches in this program
fig_width_required = 18  # cm required by the journal
cm = 2.54  # inches to cm
scaling = cm * fig_width / fig_width_required

plott(folder, collection, collection2, figsize, scaling)

dpi_target = 500  # desired dpi
dpi = dpi_target / scaling  # to get a target dpi after resizing
figname = 'Figure5.png'
path = '/nobackup/fvivant/programs/1_github/Figures/'
# path='/path/to/your/plots/'

plt.savefig('/nobackup/fvivant/programs/1_github/Figures/' + figname,
            format=figname[-3:],
            dpi=dpi,
            bbox_inches='tight')
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
