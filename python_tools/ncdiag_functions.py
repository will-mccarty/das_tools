import ncdiag as ncd
import numpy as np

def chars_to_str(chars):
    return(''.join(chars.astype(str)).strip() )

def station(data=None, return_deps=False):
    return(np.apply_along_axis(chars_to_str, axis=1, arr=data['Station_ID']))


def amb(data=None,return_deps=False):
    deps = ['omf','oma']

    val = data[ncd.var_to_var('omf')]-data[ncd.var_to_var('oma')]
    return(val)

def spd_omf(data=None,return_deps=False):
    deps = ['u_obs','v_obs','u_omf','v_omf']
    u_bkg = data[ncd.var_to_var('u_obs')] - data[ncd.var_to_var('u_omf')] 
    v_bkg = data[ncd.var_to_var('v_obs')] - data[ncd.var_to_var('v_omf')]
    spd_o = (data[ncd.var_to_var('u_obs')]**2 + data[ncd.var_to_var('v_obs')]**2) ** (0.5)
    spd_b = (u_bkg**2 + v_bkg**2) ** (0.5)

    val = spd_o - spd_b
    return(val) 

def sigo_input(data=None,return_deps=False):
    deps = ['Errinv_Input']

    val = 1.0 / data['Errinv_Input']
    msk = val > 9999.
    val[msk] = -9999.9
    return(val)

def sigo_final(data=None,return_deps=False):
    deps = ['Errinv_Final']

    val = 1.0 / data['Errinv_Final']
    msk = (val > 9999.)
    val[msk] = -9999.9
    return(val)

def qifn(data=None, return_deps=False):
    deps = ['Station_Elevation']

    qify = np.trunc(data['Station_Elevation'] / 1000.)
    qifn = data['Station_Elevation'] - (qify * 1000.)

    return(qifn)

def qify(data=None, return_deps=False):
    deps = ['Station_Elevation']

    qify = np.trunc(data['Station_Elevation'] / 1000.)
    return(qify)


def dist(data=None,return_deps=False):
    from scipy.spatial.distance import cdist

    deps = ['lon','lat']
    if (data['centroid_lat'] is None or data['centroid_lon'] is None):
       raise Exception('Distance calculations require a focal lon/lat to be set by ncd.set_centroid(lon,lat)')
    d = ( (data[ncd.var_to_var('lat')] - data['centroid_lat']) ** 2 + (data[ncd.var_to_var('lon')] - data['centroid_lon']) ** 2 ) ** (0.5) 
#    d = cdist(np.transpose( [data[ncd.var_to_var('lat'), data[ncd.var_to_var('lon')] ), [data['centroid_lat', data['centroid_lon'], 'euclidean')
    return(d)

def rand(data=None,return_deps=False):
    import numpy.random as rnd 
    deps = ['lon']

    nv = data[ncd.var_to_var('lon')].size
    val = rnd.rand(nv)
    return(val)

def cpen(val, dep=None):
    omg  = val
    sigo = dep['sigo']
    s = 0.0
    ct = 0
    for comg, csigo in zip(omg,sigo):
        ct = ct + 1
        s = s + (comg / csigo)**2
        
    return((s/ct))#**(0.5))

def haversine_pairwise(t1, t2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)

    All args must be of equal length.    

    """
    
    lon1, lat1, lon2, lat2 = map(np.radians, [t1.v('lon'), t1.v('lat'), t2.v('lon'), t2.v('lat')])

    n1 = lon1.size
    n2 = lon2.size
 
    lon1g = np.tile(lon1,(n2,1))
    lat1g = np.tile(lat1,(n2,1))
    lon2g = np.tile(np.array([lon2]).T,(1,n1))
    lat2g = np.tile(np.array([lat2]).T,(1,n1))

    lon1v = lon1g.reshape(n2*n1)
    lat1v = lat1g.reshape(n2*n1)
    lon2v = lon2g.reshape(n2*n1)
    lat2v = lat2g.reshape(n2*n1)
    
    dlon = lon2v - lon1v
    dlat = lat2v - lat1v

    a = np.sin(dlat/2.0)**2 + np.cos(lat1v) * np.cos(lat2v) * np.sin(dlon/2.0)**2

    c = 2 * np.arcsin(np.sqrt(a))
    km = 6367 * c

    ret = km.reshape(n2,n1)
    return(ret)

