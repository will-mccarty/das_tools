import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np

def plot_obs_map(cob,c=None, cmap='viridis', title=None,lat='lat',lon='lon',vmin=None,vmax=None):

#   set projection
    proj=ccrs.PlateCarree()

#   initialize the plot pointing to the projection 
    ax = plt.axes(projection=proj)

    if c is None: c='r'
#    plot the data at x=lon, y=lat.  Note I've commented out shading the dots by 
#     the air_temp value, simply because it's easier to see red
    ax.scatter(cob.v(lon),cob.v(lat),
#           c=cob.v('air_temperature@ObsValue'),
           c=c,
           linewidth=0,
           transform=proj,cmap=cmap,vmin=vmin,vmax=vmax)

#plot globally
    ax.set_global()

# draw coastlines
    ax.coastlines()

    if title is not None: 
        ax.set_title(title)

# show plot
    plt.show()

def plot_2d_hist(var1,var2,xbins=None,ybins=None,bins=40,cmin=None,cmax=None):
    if (xbins is None and ybins is not None) or (ybins is None and xbins is not None):
        raise ValueError('Currently wired for specification of both x & y bins')
    elif xbins is None and ybins is None:
        bins = bins
    else:
        bins = [xbins,ybins]

    h, xedge, yedge, im  = plt.hist2d(var1,var2,bins=bins,cmin=cmin,cmax=cmax)
    plt.show()

def plot_spatial_gauss_hist(lons,lats,deltalat=1.0,cb_label='Count',fn=None,vmin=None,vmax=None,cmap='plasma'):
    from gmao_tools import gauss_spatial_hist, gauss_spatial_hist_to_grid

    hlon, hlat, h, dlons =  gauss_spatial_hist(lons,lats,deltalat=deltalat)

    lons_new,lats_new,hnew = gauss_spatial_hist_to_grid(hlon,hlat,dlons, h, deltalat=deltalat)

    plot_spatial_gauss(lons_new,lats_new,hnew,cb_label=cb_label,fn=fn,vmin=vmin,vmax=vmax,cmap=cmap)


def plot_spatial_gauss_sum(lons,lats,val,deltalat=1.0,cb_label='Sum',fn=None,vmin=None,vmax=None,cmap='plasma'):
    from gmao_tools import gauss_spatial_sum, gauss_spatial_hist_to_grid

    hlon, hlat, h, dlons =  gauss_spatial_sum(lons,lats,val,deltalat=deltalat)

    lons_new,lats_new,hnew = gauss_spatial_hist_to_grid(hlon,hlat,dlons, h, deltalat=deltalat)

    plot_spatial_gauss(lons_new,lats_new,hnew,cb_label=cb_label,fn=fn,vmin=vmin,vmax=vmax,cmap=cmap)

def plot_spatial_gauss_mean(lons,lats,val,deltalat=1.0,cb_label='Mean',fn=None,vmin=None,vmax=None,cmap='plasma'):
    from gmao_tools import gauss_spatial_sum, gauss_spatial_hist, gauss_spatial_hist_to_grid

    hlon, hlat, s, dlons =  gauss_spatial_sum(lons,lats,val,deltalat=deltalat)
    hlon, hlat, h, dlons =  gauss_spatial_hist(lons,lats,deltalat=deltalat)

    h[h > 0] = s[h > 0]/h[h > 0]

    lons_new,lats_new,hnew = gauss_spatial_hist_to_grid(hlon,hlat,dlons, h, deltalat=deltalat)

    plot_spatial_gauss(lons_new,lats_new,hnew,cb_label=cb_label,fn=fn,vmin=vmin,vmax=vmax,cmap=cmap)


def plot_spatial_gauss(lons_new,lats_new,hnew,cb_label='',fn=None,vmin=None,vmax=None,cb_labelsize=15,cmap='plasma',suptitle=None):

    fig = plt.figure(figsize=(15,5))
    if vmin is None: vmin = np.min(hnew)
    if vmax is None: vmax = np.max(hnew)
    
    col='black'
    
    fig = plt.figure(figsize=(6.25,6))
    
    #ax0 = plt.subplot(131,projection=ccrs.Mollweide())
    ax0 = plt.subplot2grid((2, 2), (0, 0), colspan=2,projection=ccrs.Mollweide())
    ax0.pcolormesh(lons_new,lats_new,hnew,transform=ccrs.PlateCarree(),cmap=cmap,vmin=vmin,vmax=vmax)
    ax0.gridlines(linewidth=2, color=col)
    ax0.coastlines(color=col)    
    
    
    outside = (lats_new < 10)
    
    polelon = np.ma.masked_where(outside,lons_new).data
    polelat = np.ma.masked_where(outside,lats_new).data
    poleh   = np.ma.masked_where(outside,hnew).data
    
    #ax1 = plt.subplot(132,projection=ccrs.NorthPolarStereo())
    ax1 = plt.subplot2grid((2, 2), (1, 0), colspan=1,projection=ccrs.NorthPolarStereo())
    ax1.set_extent([-180, 180, 40, 90], ccrs.PlateCarree())
    ax1.pcolormesh(polelon,polelat,poleh,transform=ccrs.PlateCarree(),cmap=cmap,vmin=vmin,vmax=vmax)
    ax1.gridlines(linewidth=2, color=col)
    ax1.coastlines(linewidth=3,color=col)    
    
    
    outside = (lats_new > -10)
    
    polelon = np.ma.masked_where(outside,lons_new).data
    polelat = np.ma.masked_where(outside,lats_new).data
    poleh   = np.ma.masked_where(outside,hnew).data
    
    #ax2 = plt.subplot(133,projection=ccrs.SouthPolarStereo())
    ax2 = plt.subplot2grid((2, 2), (1, 1), colspan=1,projection=ccrs.SouthPolarStereo())
    ax2.set_extent([-180, 180, -40, -90], ccrs.PlateCarree())
    pcm = ax2.pcolormesh(polelon,polelat,poleh,transform=ccrs.PlateCarree(),cmap=cmap,vmin=vmin,vmax=vmax)
    ax2.gridlines(linewidth=2, color=col)
    ax2.coastlines(linewidth=3,color=col)    
    
    
    #fig.subplots_adjust(right=0.9)
    #cbar_ax = fig.add_axes([0.95, 0.15, 0.005, 0.7])
    #fig.colorbar(pcm,cax=cbar_ax)
    
    cbar_ax = fig.add_axes([0.865, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(pcm,cax=cbar_ax)
    cbar.ax.tick_params(labelsize=cb_labelsize)
#    cbar.set_label('Spire Occultation Profile Count\n23 Sep - 9 Dec 2018',size=13)
    cbar.set_label(cb_label,size=cb_labelsize)
    if (suptitle is not None): plt.suptitle(suptitle,fontsize=16)
    plt.tight_layout()
    
    
    plt.subplots_adjust(right=0.85)
    if fn is not None: plt.savefig(fn,bbox_inches='tight')

    plt.show()

def get_cms_vs_z(var,hgt,hgtbins=None,nbins=10):
    import gmao_tools as gt
    if hgtbins is None:
        hmin = np.min(hgt)
        hmax = np.max(hgt)

        hgtbins = gt.bins(hmin,hmax,(hmax-hmin)/nbins)

    ct = []
    mn = []
    sd = []
    sm = []
    lv = []
    lo = []
    hi = []
    for lobin,hibin in zip(hgtbins[:-1],hgtbins[1:]):
        msk = (hgt > lobin) & (hgt <= hibin)
        cct = np.sum(msk)
        cmn = np.mean(var[msk])
        csd = np.std(var[msk])
        csm = np.sum(var[msk])
        clv = (lobin+hibin)/2.

        ct.append(cct)
        mn.append(cmn)
        sd.append(csd)
        sm.append(csm)
        lv.append(clv)
        lo.append(lobin)
        hi.append(hibin)

    dct = {'ct': ct,
           'mn': mn,
           'sd': sd,
           'sm': sm,
           'lv': lv,
           'lo': lo,
           'hi': hi}

    return(dct)

def plot_cms_vs_z(ncd,var,msks,labels,colors=None,hgtbins=None,nbins=10,show=True,title=None,fn=None,scinote=[True,False,False],vdct={}):
#    vdct = {}
    ncd.use_mask(True)
    if colors is None: 
        colors = np.zeros(len(labels))# * None
    for cmsk,clab,ccol in zip(msks,labels,colors):
        ncd.set_mask(cmsk)
        vdct[clab] = get_cms_vs_z(ncd.v(var),ncd.v('Height'),hgtbins=hgtbins,nbins=nbins)
        vdct[clab]['color'] = ccol

    plot_cms_vs_z_plotter(vdct,hgtbins=hgtbins,nbins=nbins,title=title,fn=fn,scinote=scinote)
    return(vdct)

def plot_sens_vs_z(ncd,var,msks,labels,colors=None,hgtbins=None,nbins=10,show=True,title=None,fn=None,vdct={}):
    vdct = {}   
    ncd.use_mask(True)
    if colors is None:
        colors = np.zeros(len(labels))# * None
    for cmsk,clab,ccol in zip(msks,labels,colors):
        ncd.set_mask(cmsk)
        vdct[clab] = get_cms_vs_z(ncd.v(var),ncd.v('Height'),hgtbins=hgtbins,nbins=nbins)
        vdct[clab]['color'] = ccol

    plot_cms_vs_z_plotter(vdct,hgtbins=hgtbins,nbins=nbins,title=title,fn=fn,vars=['ct','mn','sm'],xlabels=['Count','FSOI per Ob.','Total FSOI'],scinote=[True,True,True])


def plot_cms_vs_z_plotter(in_vdct,hgtbins=None,nbins=10,show=True,title=None,fn=None,vars=['ct','mn','sd'],xlabels=['Count','Mean','Std. Dev'],scinote=[True,False,False]):
    import gmao_tools as gt
    import matplotlib
    from matplotlib.ticker import ScalarFormatter

    matplotlib.rcParams.update({'font.size': 16})

    fig, (axct,axmn,axsd)  = plt.subplots(1,3,sharey=True, figsize=(10.5,5))

    ctmax = 0.
    mnmax = 0.
    sdmax = 0.
    lvmax = 0.
    lvmin = 9.9e99
#    vdct = get_cms_vs_z(var,hgt,hgtbins=hgtbins,nbins=nbins)
    for cdct in in_vdct:
        print(cdct)
        print(' ')
        vdct=in_vdct[cdct]
        ct = vdct[vars[0]]
        mn = vdct[vars[1]]
        sd = vdct[vars[2]]
        lv = vdct['lv']
        lo = vdct['lo']
        hi = vdct['hi']
    
        ccol = vdct['color']  
        axct.plot(ct,lv,label=cdct,color=ccol,linewidth=2)
        axmn.plot(mn,lv,label=cdct,color=ccol,linewidth=2)
        axsd.plot(sd,lv,label=cdct,color=ccol,linewidth=2)
#        ctmax = np.max([np.nanmax(ct),ctmax])
#        mnmax = np.max([np.nanmax(np.abs(mn)),mnmax])
#        sdmax = np.max([np.nanmax(sd),sdmax])
#        print(sdmax)
        lvmax = np.max([np.nanmax(hi),lvmax])
        lvmin = np.min([np.nanmin(lo),lvmin])

        if (vars[0] == 'ct' or vars[0] == 'sd'):
            ctmax = np.max([np.nanmax(ct),ctmax])
            axct.set_xlim(0,ctmax)
        else:
            ctmax = np.max([np.nanmax(np.abs(ct)),ctmax])
            axct.set_xlim(-ctmax,ctmax)
            axct.plot([0,0],[lvmin,lvmax],'k:')
  
        if (vars[1] == 'ct' or vars[1] == 'sd'):
            mnmax = np.max([np.nanmax(mn),mnmax])
            axmn.set_xlim(0,mnmax)
        else:
            mnmax = np.max([np.nanmax(np.abs(mn)),mnmax])
            axmn.set_xlim(-1.0 * mnmax,mnmax)
            axmn.plot([0,0],[lvmin,lvmax],'k:')
        print(vars[2])
        if (vars[2] == 'ct' or vars[2] == 'sd'):
            sdmax = np.max([np.nanmax(sd),sdmax])
            axsd.set_xlim(0,sdmax)
        else:
            csdmax = np.max([np.nanmax(np.abs(sd)),sdmax])
            if (np.abs(csdmax) > sdmax): 
                sdmax = csdmax
                axsd.set_xlim(-sdmax,sdmax)
            axsd.plot([0,0],[lvmin,lvmax],'k:')

    axct.set_ylim(lvmin,lvmax)    

    offsetTextSize=12

    if (scinote[0]): 
        axct.ticklabel_format(style='sci', axis='x', scilimits=(0,0),useMathText=True)
        axct.xaxis.offsetText.set_fontsize(offsetTextSize)
    if (scinote[1]): 
        axmn.ticklabel_format(style='sci', axis='x', scilimits=(0,0),useMathText=True)
        axmn.xaxis.offsetText.set_fontsize(offsetTextSize)
    if (scinote[2]): 
        axsd.ticklabel_format(style='sci', axis='x', scilimits=(0,0),useMathText=True)
        axsd.xaxis.offsetText.set_fontsize(offsetTextSize)

    axct.xaxis.set_major_locator(plt.MaxNLocator(3))
    axmn.xaxis.set_major_locator(plt.MaxNLocator(5))
    axsd.xaxis.set_major_locator(plt.MaxNLocator(5))

    axct.set_ylabel('Height (m)')

    axct.set_xlabel(xlabels[0])
    axmn.set_xlabel(xlabels[1])
    axsd.set_xlabel(xlabels[2])
    if title is not None: plt.suptitle(title)#'Normalized Background Departure (O-F)/F')

#    axsd.legend(bbox_to_anchor=(1.675, 0.5),loc='center right', borderaxespad=0.,fontsize=14,ncol=1)

    axsd.legend(bbox_to_anchor=(1.005, 0.5),loc='center left', borderaxespad=0.,fontsize=14,ncol=1)

    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
#    axmn.xaxis.set_tick_params(rotation=90)
#    axsd.xaxis.set_tick_params(rotation=90)

    if fn is not None: plt.savefig(fn,bbox_inches='tight')
    if show: plt.show()

def plot_obs_rel_diff(stat,var,ctl,exp,in_msk='',levs=None,levs_min=None,levs_max=None,lev_var=None,color='k',minsample=0,yrange=None,fn=None,xlabel=None,ylabel=None,title=None,yscale=None,figsize=None):
    import gmao_tools as gt
    import matplotlib


    if figsize is not None: plt.figure(figsize=figsize)
    matplotlib.rcParams.update({'font.size': 16})

    if (levs is None and levs_min is None and levs_max is None):
        raise ValueError('no definition of vertical level')

    if ((levs_min is None and levs_max is not None) or 
       (levs_max is None and levs_min is not None)):
        raise ValueError('level min but not max or vice versa defined, need both')

    clevs = levs
    clevs_min = levs_min
    clevs_max = levs_max

    if (levs_max is not None and levs_min is not None and levs is None):
        clevs = (levs_max + levs_min)/2.0

    if (levs is not None and levs_min is None and levs_max is None):
        clevs_min = levs
        clevs_max = levs

    if in_msk is not '': in_msk = in_msk + ' & '

    if lev_var is None: lev_var = 'pres'
    for clev, clev_min, clev_max in zip(clevs,clevs_min,clevs_max):
        print(clev, clev_max, clev_min)
        if (levs_max is not None and levs_min is not None):
#           cmsk = inmsk + '(lev_var >= clev_min) & (lev_var < clev_max)'
            cmsk = in_msk + '({} >= {}) & ({} < {})'.format(lev_var,clev_min,lev_var,clev_max)
        else:
#           cmsk = inmsk + '(lev_var == clev)'
            cmsk = in_msk + '({} == {})'.format(lev_var,clev)

        print(cmsk)
        
        ctl.set_mask(cmsk)
        ctl.use_mask(True)
        exp.set_mask(cmsk)
        exp.use_mask(True)    

        res = gt.obs_diff_and_bars(stat,var,ctl,exp,minsample=minsample)

        plt.plot([res['rel_mean_diff']],[clev],color=color,marker='o')
        plt.plot([res['rel_mean_diff'] - res['rel_confidence_interval'],res['rel_mean_diff'] + res['rel_confidence_interval']],
                 [clev,clev],color=color,linestyle='-')

    if yrange is not None: plt.ylim(yrange)
    if ylabel is not None: plt.ylabel(ylabel)
    if xlabel is not None: plt.xlabel(xlabel)
    if title is not None:  plt.title(title)
    if yscale is not None: plt.yscale(yscale)

    plt.plot([0,0],[np.min(clevs_min),np.max(clevs_max)],'k:')
    plt.tight_layout()
    if fn is not None: plt.savefig(fn)
    plt.show()

    return
