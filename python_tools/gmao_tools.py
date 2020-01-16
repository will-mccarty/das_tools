import datetime as dt
import numpy as np

def bins(min,max,binsize):
    import numpy as np
    
    halfsize = binsize / 2.0
    binarr = np.arange(min-halfsize,max+halfsize,binsize)

    return(binarr)

def get_ndate_timeseries(sdate,edate,hr_inc=6,day_inc=0,month_inc=0,datetime=False):

       from dateutil.relativedelta import relativedelta

       sdt = ndate_to_dt(sdate)
       edt = ndate_to_dt(edate)

       tarray = []
       while (sdt <= edt):
           if (not datetime):
               tarray.append(dt_to_ndate(sdt))
           else:
               tarray.append(sdt)
#           sdt = sdt + dt.timedelta(hours=hr_inc,days=day_inc,months=mon_inc)
           sdt = sdt + relativedelta(hours=hr_inc,days=day_inc,months=month_inc)
       return(tarray)

def ndate_to_dt_single(indate):
# convert ncep ndate format (YYYYMMDDHH) to datetime object (y, m, d, h, m, s)

    yyyy   = int(indate / 1000000)
    indate = indate - (yyyy * 1000000)
    mm     = int(indate / 10000)
    indate = indate - (mm * 10000)
    dd     = int(indate / 100)
    indate = indate - (dd * 100)
    hh     = int(indate)
    return dt.datetime(yyyy, mm, dd, hh, 0, 0)

def ndate_to_dt(indate):
# convert ncep ndate format (YYYYMMDDHH) to datetime object (y, m, d, h, m, s)
    if type(indate) is int:
       out = ndate_to_dt_single(indate)
    else:

#       indate = [ indate ]

       out = []
       for date in indate:
          out.append(ndate_to_dt_single(date))

    return(out)

def dt_to_ndate(dt):
# convert datetime object to ncep ndate format (YYYYMMDDHH)
    if type(dt) is not list: 
       out = dt.year * 1000000 + dt.month * 10000 + dt.day * 100 + dt.hour
    else:
       indt = dt
       out = []
       for cdt in indt:
          out.append(cdt.year * 1000000 + cdt.month * 10000 + cdt.day * 100 + cdt.hour)

    return(out)

def geos_file_from_dt(exp,collection,dt,prefix=''):
    file = "{}.{}.{:04d}{:02d}{:02d}_{:02d}z.nc4".format(exp,collection,dt.year,dt.month,dt.day,dt.hour)
    return(prefix+file)


def find_y(x,xmin=0,xmax=0,ymin=0,ymax=0):
    y = (ymax-ymin)*(x-xmin)/(xmax-xmin)+ymin
    return(y)

def gauss_spatial_hist_to_grid(in_lons,in_lats,deltalons,hist,deltalat=1.0,out_deltall=1.0):
    import matplotlib

    lats = in_lats
    lons = in_lons
    lons[lons > 180.0] = lons[lons > 180.0] - 360.0
    
    cmap = matplotlib.cm.get_cmap('Spectral')

    nhist = (hist - np.min(hist))/(np.max(hist) - np.min(hist))
    dlat = deltalat
    nlats = lats.size
   
    nx = int(((360.0 / out_deltall) + 1))
    ny = int(((180.0 / out_deltall) + 1))
    #print(nx,ny)
    lons_new=np.arange(nx)*out_deltall-180.0
    lats_new=np.arange(ny)*out_deltall-90.0
    outlons,outlats = np.meshgrid(lons_new,lats_new)
    
    outh = np.zeros([ny,nx]) #outlons * 0.0
    
    for i,clat,clon,dlon,chist in zip(np.arange(0,nlats),lats,lons,deltalons,hist):
        cx  = (find_y(clon     ,xmin=-180,xmax=180,ymin=0,ymax=nx) ).astype(int)
        cxp = (find_y(clon+dlon,xmin=-180,xmax=180,ymin=0,ymax=nx) ).astype(int)

        cy =  (find_y(clat          ,xmin=-90,xmax=90,ymin=0,ymax=ny) ).astype(int)
        cyp = (find_y(clat+deltalat,xmin=-90,xmax=90,ymin=0,ymax=ny) ).astype(int)
        #print(cx,cxp,cy,cyp)
        outh[cy:cyp,cx:cxp] = chist
          
    return(outlons,outlats,outh)

def gauss_spatial_hist(lons,lats,deltalat=1.0):
    return(gauss_spatial_sum(lons,lats,lats*0.0+1.0,deltalat=deltalat))

def gauss_spatial_sum(in_lons,in_lats,vals,deltalat=1.0):
    histlat = np.array([])
    histlon = np.array([])
    hist    = np.array([])
    deltalons = np.array([])

    lats = in_lats
    lons = in_lons   
    lons[lons > 180.0] = lons[lons > 180.0] - 360.0
    
    for clat in np.arange(-90.,90.,deltalat):
        nlon = np.cos((clat+ clat + deltalat)/2. * np.pi / 180.0) * 360.0 / deltalat
        nlon = np.trunc(nlon)
        deltalon = 360.0 / (nlon)
        #print(nlon,deltalon)
        for clon in np.arange(-180.0,180.0,deltalon):
            msk = (lons >= clon) & (lons < clon + deltalon) & (lats >= clat) & (lats < clat + deltalat)
#            n = np.sum(msk)
            n = np.sum(vals[msk])
            histlat = np.append(histlat,clat)
            histlon = np.append(histlon,clon)
            hist = np.append(hist,n)
            deltalons= np.append(deltalons,deltalon)


    return(np.array(histlon),np.array(histlat),np.array(hist),np.array(deltalons))




def obs_diff_and_bars(stat,var,ctl,exp,minsample=0):
    import statsmodels.stats.api as sms
    import scipy.stats as st

    tsc = ctl.stat_ts(stat,var,DataFrame=True)
    tse = exp.stat_ts(stat,var,DataFrame=True)
    tsc_ct = ctl.stat_ts('count',var,DataFrame=True)
    tse_ct = exp.stat_ts('count',var,DataFrame=True)

    tsc = tsc.set_index('datetime')
    tse = tse.set_index('datetime')
    tsc_ct = tsc_ct.set_index('datetime')
    tse_ct = tse_ct.set_index('datetime')
    
    tsc = tsc.merge(tsc_ct,left_index=True,right_index=True,suffixes=('','_ct'))
    tse = tse.merge(tse_ct,left_index=True,right_index=True,suffixes=('','_ct'))
#    print(tsc)
#    tsc = tsc.mask(tsc['Val_ct'] == 0,0)
#    tse = tse.mask(tse['Val_ct'] == 0,0)
#    print(tsc)
    ts = tsc.merge(tse,left_index=True,right_index=True,suffixes=('_ctl','_exp'))

    ts = ts.mask(ts['Val_ct_ctl'] <= minsample,np.nan)
    ts = ts.mask(ts['Val_ct_exp'] <= minsample,np.nan)

    ts = ts.dropna()

    arr = (ts['Val_ctl'] - ts['Val_exp'])

    cint = st.t.interval(0.95, len(arr)-1, loc=np.mean(arr), scale=st.sem(arr))

#    cint = sms.DescrStatsW(arr,weights=ts['Val_ct_ctl']).tconfint_mean()
#    cm = sms.CompareMeans(sms.DescrStatsW(ts['Val_ctl'],weights=ts['Val_ct_ctl']/np.sum(ts['Val_ct_ctl']), sms.DescrStatsW(ts['Val_exp'],weights=ts['Val_ct_exp']/np.sum(ts['Val_ct_exp']))
#    cm = sms.CompareMeans(sms.DescrStatsW(ctl.v(var)), sms.DescrStatsW(ctl.v(var)))
#    print(ts['Val_ct_ctl'])
#    cint = cm.tconfint_diff(usevar='unequal',alpha=0.001)
#    print(cint)

    ctlmn = ts['Val_ctl'].mean()
    expmn = ts['Val_exp'].mean()
#    ctlmn = ctl.stat('mean',var)
#    expmn = exp.stat('mean',var)
#    print(ctlmn,expmn)

    delta_mean = ctlmn - expmn

    improved = False

    if (stat == 'mean'):
        if (np.abs(ctlmn) > np.abs(expmn)):
            improved = True
    else:
        if (ctlmn > expmn):
            improved = True



    dict = {
        'ctl_mean': ctlmn,
        'exp_mean': expmn,
        'abs_mean_diff': delta_mean,
        'rel_mean_diff': delta_mean / ctlmn * 100.,
        'abs_confidence_interval': (cint[1]-cint[0])/2.0,
        'rel_confidence_interval': ((cint[1]-cint[0])/2.0)/np.abs(ctlmn)*100.,
        'improved': improved
    }  

    return(dict)
