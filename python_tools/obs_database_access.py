# obs_database_access.py
#  Created: 01-2017    by:  Will McCarty (GMAO)
#
# This python module interfaces with the SEMPERPY/GMAOPY PSQL database 
# utilized on NCCS/Discover.  This module is not intended to plot the 
# observations; its aim is to directly interface/return values and statistics
# from the database so that the user can generate their own plots (or wrappers)
#
# Dependencies:  Numpy, copy, pprint, psycopg2 
#

import numpy as np
import datetime as dt

def query_structure(date=None,variable=None,domain_name=None,level=None,levtype=None,kt=None,kx=None, usage=None):

# function query_structure()
#
#  Description:  Creates a dictionary that can be used for querying the 
#    database.  This query is passed to the functions contained within the 
#    experiment class below.  
#  Notes:
#       - The list of 'expected' is not exhaustive.  See 
#           experiment.show_all_distinct() below.  The inputs represent 
#           initial values that can be set in the returned dictionary.  
#       - Once created, the dictionary can be modified by hand (see example).
#       - The 'statistic' DB table field is explicitly excluded here, as it is 
#           wired relative to the statistic query/generation routines below.  
#           If there were ever a reason to include it, it can be added as its 
#           own key by hand (in theory - never tested)
#       - Date handling is intended to be fairly special.  If dates are 
#           explicitly specified in the query, then those dates will be 
#           used.  If not, each experiment instance (see the class below)
#           will have a start and end date set.  It will assume that all
#           values to be queried are to be within the start/end date if the
#           date key in the query is None/does not exist.
#
#  Inputs:       (default     |   expected)        
#    date:        None            YYYYMMDDHH (individual or list)
#    variable:    None            'omf', 'oma'
#    domain_name: None            'n.hem', 's.hem', 'global'
#    level:       None            200., 300., 850.
#    levtype:     None            'pres', 'sfc'
#    kt:          None            4, 5, 11 (see gmao KT list)
#    kx:          None            120, 180, 220 (see gmao KX list)
#    usage:       None            '37', 'used', 'passive
#  Returns:
#    Dictionary containing the above entries.  
#  Example:
#    import obs_database_access as obsdb
#    from pprint import pprint
#
#    query = obsdb.query_structure(date=[2006010100, 2006010106], 
#                                  variable='omf')
#    pprint(query)
#
#    query['variable'] = 'oma'
#    pprint(query)

    query = {
      'date':        date        ,
      'variable':    variable    ,
      'domain_name': domain_name ,
      'level':       level       ,
      'levtype':     levtype     ,
      'kx':          kx          ,
      'kt':          kt          ,
      'usage':       usage       }
    return(query)

##def ndate_to_dt(indate):
## convert ncep ndate format (YYYYMMDDHH) to datetime object (y, m, d, h, m, s)
#    yyyy   = int(indate / 1000000)
#    indate = indate - (yyyy * 1000000)
#    mm     = int(indate / 10000)
#    indate = indate - (mm * 10000)
#    dd     = int(indate / 100)
#    indate = indate - (dd * 100)
#    hh     = int(indate)
#    return dt.datetime(yyyy, mm, dd, hh, 0, 0)
#
#def dt_to_ndate(dt):
## convert datetime object to ncep ndate format (YYYYMMDDHH)
#    ndate = dt.year * 1000000 + dt.month * 10000 + dt.day * 100 + dt.hour
#
#    return ndate


class experiment():
# Class experiment 
#
#  Description:  Class (interface and subroutines) used to connect to and query
#    the SEMPERPY/GMAOPY database on discover. 
#  Notes: 
#    - Currently, the default is the user database, though it should work with 
#      the OPS database as long as it is explicitly specified in the 
#      initialization - optional flag database='ob_ops' (I think).
#  Inputs (upon initialization):
#    exp:       experiment ID used to populate the database
#    database:  database name (optional, ob_exp by default)
#    startdate: YYYYMMDDHH start date (optional, will determine first date for
#                  experiment otherwise)
#    enddate:   YYYYMMDDHH end date (optional, will determine last date for
#                  experiment otherwise)
#  Returns:
#    self instance: variable containing the conenction and used to query
#                    for initialized experiment
#  Usage:
#    import obs_database_access as obsdb
#    from pprint import pprint
#
#    query = obsdb.query_structure(variable='omf', domain='global,
#                                  usage='used', level=200., kt=[4,5]
#                                  kx=220)
#    exp = obsdb.experiment('m2m_1pch',startdate=2006070100,enddate=2006070218)
#    
#    #Now you are connected and ready to query the DB via variable exp.  
#    #This variable/connection will be used in the examples below
    def __init__(self, exp, database='ob_exp', startdate=None, enddate=None, verbose_query=None):
        import psycopg2 as psql

        self.experiment = exp
        print('Initializing as ',exp)
        self.database = database
        if (database == 'ob_exp'):
           self.connect_string = "dbname='semper' user='gmao_user' host='edb1' "
        elif (database == 'ob_ops'):
           self.connect_string = "dbname='gmao_stats' user='gmao_user' host='edb1' "
        else:
           print('Warning, not DB ob_exp or ob_ops, using dbname \'semper\'')
           self.connect_string = "dbname='semper' user='gmao_user' host='edb1' "

        self.con = psql.connect(self.connect_string)
        if (startdate and enddate):
            self.startdate = startdate
            self.enddate = enddate
        elif (startdate or enddate):
            print('warning, only startdate or enddate set, determining it from DB')
            self.startdate, self.enddate = self.get_min_max_date()
        else:
            self.startdate, self.enddate = self.get_min_max_date()
        print('experiment start, end dates:',self.startdate, self.enddate)
        self.fields = ['north' ,'west'  ,'south'   ,'east'     ,'domain_name','level','levtype',
                      'expver','source','variable','statistic','type'       ,'kx'   ,'kt'     ,
                      'usage' ,'date'  ,'value'   ,'count']
        self.verbose_query = verbose_query

    def set_daterange(self, startdate, enddate):
# Function exp.set_daterange()
#  Description:  Alter the start/end dates set upon initialization.
#  Notes:
#    - Should complain if only one date is set.  Perhaps in the future
#      make dates optional so that only one can be adjusted on the fly
#  Inputs: 
#    startdate: YYYYMMDDHH start date 
#    enddate:   YYYYMMDDHH end date 
#  Returns: 
#    Nothing
#  Example:
#    #continuing from initialization example above)
#    exp.set_daterange(2006070100,2006073118)
        self.startdate = startdate
        self.enddate = enddate
        print('experiment start, end dates set to:',self.startdate, self.enddate)

    def get_min_max_date(self,checkall=False,checkkx=220,checklev=1000,checkkt=4):
# Function exp.get_min_max_date()()
#  Description:  Determines the earliest and latest date for the experiment
#    in the database.
#  Notes:
#    - Used internally to set start/end date if not specified upon 
#        initialization
#  Inputs: 
#    Nothing
#  Returns: 
#    startdate:  YYYYMMDDHH for earliest date
#    enddate:    YYYYMMDDHH for latest date
#  Example:
#    #continuing from initialization example above)
#    start_date, end_date = exp.get_min_max_date() 
	    
        cur = self.con.cursor()
        if (checkall):
            query = "SELECT MIN(date), MAX(date) from {}.v_view where  expver = \'{}\'".format(self.database, self.experiment)
        else:
            query = "SELECT MIN(date), MAX(date) from {}.v_view where  expver = \'{}\' and kt = \'{}\' and kx = \'{}\' and level = \'{}\'".format(self.database, self.experiment,checkkt,checkkx,checklev)

        cur.execute(query)
        val = cur.fetchall()
        val = val[0]
        startdate, enddate = val
        return(startdate, enddate)

    def append_query_daterange(self):
# Function exp.append_query_daterange()
#  Description:  Returns a string to be appended to a SQL query limiting the 
#    query to be >= internally set startdate and <= internally set enddate
#  Notes:
#    - Used internally to limit date range if no date is set in a query dict
#  Inputs: 
#    Nothing
#  Returns:
#    str_to_append:  string 'and date >= startYYYYMMDDHH and 
#                            date <= endYYYYMMDDHH)'
#  Example:
#    #continuing from initialization example above)
#    query_str = 'dummy wont work'
#    query_str = query_str + exp.append_query_daterange()
#    print(query_str)
	    
        str = 'and date >= {} and date <= {} '.format(self.startdate, self.enddate)      
        return(str)

    def append_query(self, query):
# Function exp.append_query()
#  Description:  Returns a string to be appended to a SQL query based on query 
#    fields set in query structure
#  Notes:
#    - Used internally to append SQL queries with 'where' statements 
#        corresponding to the query dict
#    - Uses in statements for lists in query dict, otherwise explicit equals
#        if only one variable in a key is to be requested
#  Inputs: 
#    query:  query dict
#  Returns:
#    str_to_append:  string 'where key = val and key2 in [array,of,vals]'
#  Example:
#    #continuing from initialization example above)
#    query_str = 'dummy wont work'
#    query_str = query_str + exp.append_query(query)
#    print(query_str)
        str = ''
        hasDate = False
        for key in query:
            if (query[key]):
                if key == 'date':
                    hasDate = True
                if isinstance(query[key], list):
                    instr = 'and {} in ('.format(key)
                    init=True
                    for q in query[key]:
                        if (init):
                            init=False
                        else:
                            instr = instr + ','
                        instr = instr + '{}'.format(q)
   
                    str = str + instr + ') '
                else:
                    str = str + " and {} = \'{}\' ".format(key,query[key])
        if not hasDate:
            str = str + self.append_query_daterange()

        return(str)

    def get_values_counts(self, querydict, statistic=None):
# Function exp.get_values_counts()
#  Description:  The fundemental query of the SEMPERPY/GMAOPY database.  All 
#    metrics are saved generally as a value and a count.  e.g., for RMS, the
#    statistics is 'sumsquare', and to determine RMS, you would add all the 
#    values from a query, sum(values), divide by the total count, sum(counts)
#    and take the root.  
#  Notes:
#    - So far there is no routine to only get counts
#    - As of documentation, I've only done mean and RMS
#  Inputs: 
#    query:  query dict to retrieve values and counts arrays
#    statistic:  optional string to get only a specific statistic 
#        (e.g. sum, sumsquare)
#  Returns:
#    values:  numpy array of values for all fields matching the query
#    counts:  numpy array of the counts corresponding to the values retrieved
#  Example:
#    #continuing from initialization example above)
#    val, ct = query_str + exp.get_values_counts(query,statistic='sumsquare')
#    import numpy as np
#    print(np.sum(ct))
#    # the total number of observations in the database that fit the query
#    # will be printed
        cur = self.con.cursor()
        query = "SELECT value, count from {}.v_view where expver = \'{}\' ".format(self.database, self.experiment)
        query = query + self.append_query(querydict)
        if statistic:
            query = query + " and statistic = \'{}\' ".format(statistic)

        if (self.verbose_query is not None):
            print(query)

        cur.execute(query)
        out = cur.fetchall()
        val = []
        ct  = []

        for o in out:
            val.append(o[0])
            ct.append(o[1]) 
         
        return(val, ct)      

    def get_count(self,query):
        val, ct = self.get_values_counts(query, statistic='sum')

        count = np.sum(ct)

        return(count)

    def get_rms(self,query):
# Function exp.get_rms()
#  Description:  Get the RMS from the database for those observations that match
#    the inputted query dict
#  Notes:
#    - It is a good idea to run the show_all_distinct on a query to make sure 
#        you are really querying exactly what intend.  It will show you all
#        that will be selected, and it will prevent you from mistakenly 
#        including things
#  Inputs: 
#    query:  query dict
#  Returns:
#    rms:  float of RMS
#  Example:
#    #continuing from initialization example above)
#    print(exp.get_rms(query))
        val, ct = self.get_values_counts(query, statistic='sumsquare')

        rms = (np.sum(val) / np.sum(ct))**(0.5)
       
        return(rms)

    def get_bias(self,query):
# Function exp.get_bias()
#  Description:  Get the bias from the database for those observations that
#    match the inputted query dict
#  Notes:
#    - It is a good idea to run the show_all_distinct on a query to make sure 
#        you are really querying exactly what intend.  It will show you all
#        that will be selected, and it will prevent you from mistakenly 
#        including things
#  Inputs: 
#    query:  query dict
#  Returns:
#    bias:  float of bias
#  Example:
#    #continuing from initialization example above)
#    print(exp.get_bias(query))
        val, ct = self.get_values_counts(query, statistic='sum')

        bias = np.sum(val) / np.sum(ct)

        return(bias)

    def get_rms_from_list(self, query, field, list):
# Function exp.get_rms_from_list()
#  Description:  Get the rms from the database for those observations that
#    match the inputted query dict, except for the field inputted.  That field
#    will be looped over the list values also inputted
#  Notes:
#    - It is a good idea to run the show_all_distinct on a query to make sure 
#        you are really querying exactly what intend.  It will show you all
#        that will be selected, and it will prevent you from mistakenly 
#        including things
#  Inputs: 
#    query:  query dict
#    field:  string of the key/field in query that you want to alter and loop
#              over
#    list:   values over which you want to loop
#  Returns:
#    rms_array:  numpy array of rms floats of size (list)
#  Example:
#    #continuing from initialization example above)
#    levs_to_compute = [200.,500.,850.]
#    rms_arr = exp.get_rms_from_list(query,level,lev_to_compute)
#    for lev, rms in zip(levs_to_compute,rms_arr):
#       print(lev,rms)
#
        from copy import deepcopy

        in_query = deepcopy(query)

        rms_array = []
        for l in list:
            in_query[field] = l
            rms = self.get_rms(in_query)
            rms_array.append(rms)

        return(np.array(rms_array))

    def get_bias_from_list(self, query, field, list):
# Function exp.get_bias_from_list()
#  Description:  Get the bias from the database for those observations that
#    match the inputted query dict, except for the field inputted.  That field
#    will be looped over the list values also inputted
#  Notes:
#    - It is a good idea to run the show_all_distinct on a query to make sure 
#        you are really querying exactly what intend.  It will show you all
#        that will be selected, and it will prevent you from mistakenly 
#        including things
#  Inputs: 
#    query:  query dict
#    field:  string of the key/field in query that you want to alter and loop
#              over
#    list:   values over which you want to loop
#  Returns:
#    bias_array:  numpy array of bias floats of size (list)
#  Example:
#    #continuing from initialization example above)
#    levs_to_compute = [200.,500.,850.]
#    bias_arr = exp.get_bias_from_list(query,level,lev_to_compute)
#    for lev, bias in zip(levs_to_compute,bias_arr):
#       print(lev,bias)
#
        from copy import deepcopy
 
        in_query = deepcopy(query)

        bias_array = []
        for l in list:
            in_query[field] = l
            bias = self.get_bias(in_query)
            bias_array.append(bias)

        return(np.array(bias_array))

    def get_count_from_list(self, query, field, list):
        from copy import deepcopy

        in_query = deepcopy(query)

        count_array = []
        for l in list:
            in_query[field] = l
            count = self.get_count(in_query)
            count_array.append(count)

        return(np.array(count_array))


    def show_all_distinct(self,query):
# Function exp.show_all_distinct()
#  Description:  Loop over a query dict and show all distinct values for each
#    field/key in the database that match the query.
#  Notes:
#    - This is very useful to show what may accidentally be included
#    - This may take a long time if you provide a blankly initialized
#        query dict as it will tell you every value for every key
#    - Right now it does values and count, though that may be made optional,
#        since the values are generally all distinct
#  Inputs: 
#    query:  query dict
#  Returns:
#    Nothing (but prints output)
#  Example:
#    #continuing from initialization example above)
#    exp.show_all_distinct(query)
        import pprint

        print('Showing all distinct values matching query:')     
        pprint.pprint(query)
        
        for f in self.fields:
            d = self.get_distinct(f, in_query=query)
            print(f,d)
 
    def get_distinct(self, field, in_query=None):
# Function exp.get_distinct()
#  Description:  Returns an array of all distinct values for a specified field
#    and query, optionally.
#  Notes:
#    - If no query is provided and the date range internally is set long, this
#        may take a long time
#  Inputs:
#    field:  field/key to get distict/unique values of
#    query:  optional query dict to limit search for distinct values
#  Returns:
#    arr_distinct:  array of distinct values from the database
#  Example:
#    #continuing from initialization example above)
#    newquery = obsdb.query_structure(variable='omf', domain='global,
#                                  usage='used', kt=[4,5],kx=220)
#    distinct_levels = exp.get_distinct('level',in_query=newquery)
#    print(query_str)
        cur = self.con.cursor()
        query = "SELECT DISTINCT {} from {}.v_view where expver = \'{}\' ".format(field, self.database, self.experiment) 
#        query = "SELECT DISTINCT {} from {}.v_view where expver = \'{}\' ".format(field, 'ob_exp', self.experiment) 
       
        if (in_query):
            query = query + self.append_query(in_query)
        else:
            query = query + self.append_query_daterange()

        cur.execute(query)
        val = cur.fetchall()
        
        return(val)  

    def get_ndate_timeseries(self,sdate=None,edate=None,inc=6,datetime=False):
       import gmao_tools as gt

       if (sdate is None): sdate = self.startdate
       if (edate is None): edate = self.enddate
     
       tarray = gt.get_ndate_timeseries(sdate,edate,hr_inc=inc,datetime=datetime)
#       sdt = ndate_to_dt(sdate)
#       edt = ndate_to_dt(edate)
#
#       tarray = []
#       while (sdt <= edt):
#           if (not datetime):
#               tarray.append(dt_to_ndate(sdt))
#           else:
#               tarray.append(sdt)
#           sdt = sdt + dt.timedelta(hours=inc)

       return(tarray)
       


