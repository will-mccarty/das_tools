import datetime as dt

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
    if type(indate) is int: indate = [ indate ]
    out = []
    for date in indate:
        out.append(ndate_to_dt_single(date))

    return(out)

def dt_to_ndate(dt):
# convert datetime object to ncep ndate format (YYYYMMDDHH)
    ndate = dt.year * 1000000 + dt.month * 10000 + dt.day * 100 + dt.hour

    return ndate

