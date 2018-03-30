import ncdiag as ncd

def amb(data=None,return_deps=False):
    deps = ['omf','oma']

    val = data[ncd.var_to_var('omf')]-data[ncd.var_to_var('oma')]
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

def cpen(val, dep=None):
    omg  = val
    sigo = dep['sigo']
    s = 0.0
    ct = 0
    for comg, csigo in zip(omg,sigo):
        ct = ct + 1
        s = s + (comg / csigo)**2
        
    return((s/ct))#**(0.5))
