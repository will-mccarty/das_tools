import datetime as dt

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
