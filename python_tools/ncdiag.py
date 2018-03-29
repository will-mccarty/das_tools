import numpy as np
import datetime as dt
import netCDF4 as nc4
import ncdiag_functions as ncf

varmap = { 
    'lat':   'Latitude'         ,
    'lon':   'Longitude'        ,
    'kx':    'Observation_Type' ,
    'used':  'Analysis_Use_Flag',
    'pres':  'Pressure',
    'ob':    'Observation',
    'obs':   'Observation',
    'omf':   'Obs_Minus_Forecast_adjusted',
    'omfbc': 'Obs_Minus_Forecast_adjusted',
    'omfnbc':'Obs_Minus_Forecast_unadjusted',
    'oma':   'Obs_Minus_Analysis_adjusted',
    'omabc': 'Obs_Minus_Analysis_adjusted',
    'omanbc':'Obs_Minus_Analysis_unadjusted',
    'ichan': 'Channel_Index',
    'qcmark': 'QC_Flag',
    'chused': 'use_flag'}

derived_var = {
    'amb':         {'func': ncf.amb,          'deps': ['omf','oma']},
    'sigo_input':  {'func': ncf.sigo_input,   'deps': ['Errinv_Input']},
    'sigo_final':  {'func': ncf.sigo_final,   'deps': ['Errinv_Final']},
    'sigo':        {'func': ncf.sigo_final,   'deps': ['Errinv_Final']} 
    }

def var_to_var(in_var):
    if (in_var in varmap):
        var = varmap[in_var]
    else:
        var = in_var
    return(var)

def getFields(text):
    import re

    wds = re.compile('\w+\.*\w*').findall(text)
    return(wds[::2])


class obs():

    def __init__(self, fn, date=None, mask=None, verbose=False, reallyverbose=False):
        self.fn = fn
        self.nc4 = nc4.Dataset(self.fn)
        self.data = { 'fn':   self.fn,
                      'date': date   } 
        if mask is None:
            self.mask = None
        else:
            self.set_mask(mask)

        self.mask_enabled = None
        self.verbose = verbose
        self.reallyverbose = reallyverbose


    def use_mask(self, msk):
        self.mask_enabled = msk

    def set_mask(self, logic):
        flds = getFields(logic)

        newlogic = logic
        if (self.reallyverbose): print(newlogic)

        for fld in list(set(flds)):
            cur = self.v(fld)
            newlogic = newlogic.replace(fld,'self.data[\'{}\']'.format(var_to_var(fld)))

        if (self.reallyverbose): print(newlogic)

        self.mask=eval(newlogic)

    def v(self, in_var,masked=None):
        if masked is True and self.mask is None:
            raise ValueError('Masked is asked for, but no mask is set - use self.set_mask')

        if (masked is None): 
            msk = self.mask_enabled
        else:
            msk = masked

        derived = in_var in derived_var
        if (derived):
            vars = derived_var[in_var]['deps']
        else:
            vars = [in_var]

        for cvar in vars:
            var = var_to_var(cvar)
            try:
                self.data[var] = self.nc4.variables[var][...]
            except:
                raise ValueError('Field {} not in file'.format(var))

        var = var_to_var(in_var)
        if (derived):
            self.data[var] = derived_var[var]['func'](data=self.data)

        if (msk):
            return(self.data[var][self.mask])
        else:
            return(self.data[var])

class obs_template():
    

    def __init__(self, fn_tmpl, startdate=None, enddate=None, verbose=False, reallyverbose=False):
        import gmao_tools as gt

        self.fn_tmpl   = fn_tmpl
        if startdate is not None: self.startdate = startdate
        if enddate is not None:   self.enddate   = enddate
        self.verbose = verbose
        self.reallyverbose = reallyverbose 

        self.mask_logic = None
        self.mask_enabled = None

        self.obdict = {}

    def set_mask(self, logic):
        self.mask_logic = logic

    def use_mask(self, msk):
        self.mask_enabled = msk

    def fn_from_tmpl(self, dattim):
        from string import Template

        fn = Template(self.fn_tmpl)
        res = fn.safe_substitute(yyyy=dattim.strftime('%Y'), 
                        mm=dattim.strftime('%m'),
                        dd=dattim.strftime('%d'),
                        hh=dattim.strftime('%H'))
        return(res)


    def v(self, in_var, startdate=None, enddate=None, dates=None, hr_inc=6, masked=None):
        import gmao_tools as gt

        if (masked or self.mask_enabled) and self.mask_logic is None:
            raise ValueError('Masked data requested, but mask logic is not set')
        elif ((masked or self.mask_enabled) and self.mask_logic is not None):
            msk = self.mask_logic
        else:
            msk = None

        if startdate is None: startdate=self.startdate
        if enddate is None:   enddate=self.enddate 

        if (startdate is None and enddate is None and dates is None):
            raise ValueError('No dates are specified, unable to retrieve data')
        elif (startdate is not None and enddate is not None and dates is not None):
            raise ValueError('Start/end dates and dates array specified, unable to retrieve data due to conflicting date specifications')
        elif ( (startdate is not None and enddate is None) or (startdate is None and enddate is not None) ):
            raise ValueError('Start/end date specified without the other')
        elif (startdate is not None and enddate is not None and dates is None):
            ndates = gt.get_ndate_timeseries(startdate,enddate,hr_inc=hr_inc)
        elif ( (startdate is None and enddate is None) and dates is not None):
            ndates = dates
        else:
            raise ValueError('Unexcepted failure on date inputs, startdate = {}, enddate = {}, dates = {}'.format(startdate,enddate,dates))

        dts = gt.ndate_to_dt(ndates)

        data = None

        for dt in dts:
            fn = self.fn_from_tmpl(dt)

            if fn not in self.obdict: 
                if self.verbose: print('Opening {}'.format(fn))
                self.obdict[fn] = obs(fn,date=dt,verbose=self.verbose, reallyverbose=self.reallyverbose)
            if msk is not None:
                self.obdict[fn].set_mask(msk)

            if (data is None):
                data = self.obdict[fn].v(in_var,masked=msk)
            else:
                data = np.append(data,self.obdict[fn].v(in_var,masked=msk),axis=0)
                
        return(data)                  
