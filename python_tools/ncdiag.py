import numpy as np
import datetime as dt
import netCDF4 as nc4
import ncdiag_functions as ncf

varmap = { 
    'lat':   'Latitude'         ,
    'lon':   'Longitude'        ,
    'kx':    'Observation_Type' ,
#    'used':  'Analysis_Use_Flag',
    'pres':  'Pressure',
    'ob':    'Observation',
    'obs':   'Observation',
    'omf':   'Obs_Minus_Forecast_adjusted',
    'omfbc': 'Obs_Minus_Forecast_adjusted',
    'omfnbc':'Obs_Minus_Forecast_unadjusted',
    'oma':   'Obs_Minus_Analysis_adjusted',
    'omabc': 'Obs_Minus_Analysis_adjusted',
    'omanbc':'Obs_Minus_Analysis_unadjusted',
    'u_obs':  'u_Observation',
    'u_omf':   'u_Obs_Minus_Forecast_adjusted',
    'u_omfbc': 'u_Obs_Minus_Forecast_adjusted',
    'u_omfnbc':'u_Obs_Minus_Forecast_unadjusted',
    'v_obs':  'v_Observation',
    'v_omf':   'v_Obs_Minus_Forecast_adjusted',
    'v_omfbc': 'v_Obs_Minus_Forecast_adjusted',
    'v_omfnbc':'v_Obs_Minus_Forecast_unadjusted',
    'ichan': 'Channel_Index',
    'qcmark': 'QC_Flag',
    'chused': 'use_flag',
    'ee':     'Height'}    # Note:  amv expected error is stored in Height

derived_var = {
    'amb':         {'func': ncf.amb,          'deps': ['omf','oma']},
    'sigo_input':  {'func': ncf.sigo_input,   'deps': ['Errinv_Input']},
    'sigo_final':  {'func': ncf.sigo_final,   'deps': ['Errinv_Final']},
    'sigo':        {'func': ncf.sigo,         'deps': ['Errinv_Final','Inverse_Observation_Error']},
    'dist':        {'func': ncf.dist,         'deps': ['lon','lat']},
    'rand':        {'func': ncf.rand,         'deps': ['lon']},
    'station':     {'func': ncf.station,      'deps': ['Station_ID']},
    'spd_omf':     {'func': ncf.spd_omf,      'deps': ['u_obs','v_obs','u_omf','v_omf']},
    'qifn':        {'func': ncf.qifn,         'deps': ['Station_Elevation']},
    'qify':        {'func': ncf.qify,         'deps': ['Station_Elevation']},
    'used':        {'func': ncf.used,         'deps': ['use_flag','Analysis_Use_Flag','QC_Flag','Channel_Index']}
    }

stats = {
    'mean':        {'func': np.mean,          'deps': None} ,
    'cpen':        {'func': ncf.cpen,         'deps': ['sigo']}
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

    def __init__(self, fn, date=None, mask=None, verbose=False, reallyverbose=False, in_data=None):
        self.fn = fn
        try:
            self.nc4 = nc4.Dataset(self.fn)
        except:
            raise ValueError('NCDIAG File {} failed'.format(self.fn))

        self.data = { 'fn':   self.fn,
                      'date': date   ,
                      'centroid_lat': None,
                      'centroid_lon': None}
        for key in derived_var:
            self.data[key] = None

        self.centroid_calculated = False
        

        if in_data is not None:
           for key in in_data:
               self.data[key] = in_data[key]

        if mask is None:
            self.mask = None
        else:
            self.set_mask(mask)

        self.mask_enabled = None
        self.verbose = verbose
        self.reallyverbose = reallyverbose

    def set_centroid(self, lon, lat):
        self.data['centroid_lon'] = lon
        self.data['centroid_lat'] = lat 
        self.data['dist'] = None

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
        import warnings as warn

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
#            try:
            if (var not in self.data and var in self.nc4.variables):
                self.data[var] = self.nc4.variables[var][...]
#            except:
#                raise ValueError('Field {} not in file'.format(var))

        var = var_to_var(in_var)
        
        if (derived and self.data[var] is None):
            self.data[var] = derived_var[var]['func'](data=self.data)

        if (msk):
            return(self.data[var][self.mask])
        else:
            return(self.data[var])

    def stat(self,stat,var,masked=None):
        if masked is True and self.mask is None:
            raise ValueError('Masked is asked for, but no mask is set - use self.set_mask')

        if (masked is None):
            msk = self.mask_enabled
        else:
            msk = masked
 
        dep_dict = {}

        if stats[stat]['deps'] is not None:
            for dep in stats[stat]['deps']:
                dep_dict[dep] = self.v(dep)
            value = stats[stat]['func'](self.v(var),dep=dep_dict)
        else:
            value = stats[stat]['func'](self.v(var)) 

        return(value)

class obs_template(obs):
    

    def __init__(self, fn_tmpl, startdate=None, enddate=None, verbose=False, reallyverbose=False, hr_inc=6):
        import gmao_tools as gt

        self.fn_tmpl   = fn_tmpl
        if startdate is not None: self.startdate = startdate
        if enddate is not None:   self.enddate   = enddate
        self.hr_inc=hr_inc
        self.verbose = verbose
        self.reallyverbose = reallyverbose 

        
        self.mask_logic = None
        self.mask_enabled = None

        self.obdict = {}

        self.data = {}

    def set_mask(self, logic):
        self.mask_logic = logic

# use_mask inhereted from obs

    def fn_from_tmpl(self, dattim):
        from string import Template

        fn = Template(self.fn_tmpl)
        res = fn.safe_substitute(yyyy=dattim.strftime('%Y'), 
                        mm=dattim.strftime('%m'),
                        dd=dattim.strftime('%d'),
                        hh=dattim.strftime('%H'))
        return(res)

    def set_centroid(self, lon, lat):
        self.data['centroid_lon'] = lon
        self.data['centroid_lat'] = lat
        for cobd in self.obdict:
            self.obdict[cobd].set_centroid(lon,lat)

    def v(self, in_var, startdate=None, enddate=None, dates=None, hr_inc=None, masked=None):
        import gmao_tools as gt

        if (masked or self.mask_enabled) and self.mask_logic is None:
            raise ValueError('Masked data requested, but mask logic is not set')
        elif ((masked or self.mask_enabled) and self.mask_logic is not None):
            msk = self.mask_logic
        else:
            msk = None

        

        if startdate is None: startdate=self.startdate
        if enddate is None:   enddate=self.enddate 
        if hr_inc is None:    hr_inc=self.hr_inc

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
                self.obdict[fn] = obs(fn,date=dt,verbose=self.verbose, reallyverbose=self.reallyverbose, in_data=self.data)
            if msk is not None:
                self.obdict[fn].set_mask(msk)

            if (data is None):
                data = self.obdict[fn].v(in_var,masked=msk)
            else:
                data = np.append(data,self.obdict[fn].v(in_var,masked=msk),axis=0)
                
        return(data)                  
