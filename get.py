import glob
import icecap as icp
import numpy as np
import os
import fnmatch
import subradar as sr
import rsr
import pandas as pd


def params():
    """get various parameters defining the season
    """
    out = {'code_path':os.getcwd()}
    out['season'] = out['code_path'].split('/')[-3]
    out['process'] = out['code_path'].split('/')[-1]
    out['root_path'] = '/'.join(out['code_path'].split('/')[0:-5])
    out['norm_path'] = out['root_path'] + '/targ/norm'
    out['rsr_path'] = out['code_path'].replace('code', 'targ')
    out['cmp_path'] = out['rsr_path'].replace('RSR', 'CMP')
    out['pik_path'] = out['root_path'] + '/orig/xtra/'+out['season']+'/PIK/' + out['process']
    out['foc_path'] = out['root_path'] + '/targ/xtra/' + out['season']+ '/FOC/Best_Versions/S1_POS'
    out['sweep_path'] = out['root_path'] + '/targ/xtra/' + out['season']+ '/FOC/Best_Versions/S5_VEW'
    out['tpro_path'] = out['root_path'] + '/targ/tpro'
    out['treg_path'] = out['root_path'] + '/targ/treg'
    out['season_flight_pst'] = out['root_path'] + '/syst/linux/lib/dbase/season_flight_pst'
    return out


def pik(pst, process=None, **kwargs):
    """Get available PIK files for a PST
    """
    p = icp.get.params()
    if process is None:
        process = p['process']
    folder = '/'.join([p['pik_path'].replace('/'+p['process'],''), process, pst])
    files = glob.glob(folder + '/*.*')
    names = [i.split('/')[-1] for i in files]
    products = [i.split('.')[0] for i in names]
    pik = products
    #pik = [i.split('.')[1] for i in names]
    return products, pik


def cmp(pst, process=None, **kwargs):
    """Get available radar data in CMP for a PST
    """
    p = icp.get.params()
    if process is None:
        process = p['process']
    folder = '/'.join([p['cmp_path'].replace('/'+p['process'],''), process, pst])
    files = glob.glob(folder + '/*[!.meta]')
    products = [i.split('/')[-1] for i in files]
    return products


def pst(pattern, **kwargs):
    """Get PSTs for the current season that match a given pattern (regex)
    """
    p = icp.get.params()
    data = np.genfromtxt(p['season_flight_pst'], delimiter=' ', dtype=np.str)
    i = np.where(data[:,2] == p['season'])
    pst = data[i,0]
    return fnmatch.filter(pst.flatten(), pattern)


def sweep(pst, **kwargs):
    """Get available sweeps files for a PST
    """
    p = icp.get.params()
    folder = '/'.join([p['sweep_path'], pst])
    files = glob.glob(folder + '/*sweeps*')
    products = [i.split('/')[-1] for i in files]
    return products


def rsr(pst, process=None, **kwargs):
    """Get available rsr files
    """
    p = icp.get.params()
    if process is None:
        process = p['process']
    folder = '/'.join([p['rsr_path'].replace('/'+p['process'],''), process, pst])
    files = glob.glob(folder + '/*.*')
    products = [i.split('/')[-1] for i in files]
    pik = [i.split('.')[1] for i in products if len(i.split('.')) == 2]
    return pik


def rsr_data(pst, **kwargs):
    """Display data avaialble to launch RSR
    """
    psts = icp.get.pst(pst)
    cmps = [ icp.get.cmp(i, process='pik1') for i in psts ]
    cmps_1m = [ icp.get.cmp(i, process='pik1.1m') for i in psts ]
    piks = [ icp.get.pik(i, process='pik1')[1] for i in psts]
    piks_1m = [ icp.get.pik(i, process='pik1.1m')[1] for i in psts]
    sweeps = [ icp.get.sweep(i) for i in psts ]
    rsr_1m = [ icp.get.rsr(i, process='pik1.1m') for i in psts  ]
    d = {'PST':psts}
    df = pd.DataFrame(d)
    #df['CMP_pik1'] = cmps
    df['sweeps'] = sweeps
    df['CMP_pik1.1m'] = cmps_1m
    df['PIK_pik1'] = piks
    df['PIK1_pik1.1m'] = piks_1m
    df['RSR_1m'] = rsr_1m
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)
    return df


def flight(pst):
    """Get Flight for a PST
    """
    p = icp.get.params()
    data = np.genfromtxt(p['season_flight_pst'], delimiter=' ', dtype=np.str)
    i = np.where(data[:,0] == pst)
    return data[i,:].flatten()[1]


def surface_range(pst, **kwargs):
    """Range to surface interpolated alon the foc_time
    """
    p = icp.get.params()
    t, val = icp.read.norm(pst, 'LAS', 'las_rng', interp=True)
    return val


def ice_thickness(pst, **kwargs):
    p = icp.get.params()
    tref = icp.read.ztim(p['foc_path']+'/'+pst+'/ztim_DNhH')['htim']
    try:
        t, val = icp.read.targ(pst, 'treg', 'TRJ_JKB0', 'ztim_llzrphsaaa', interp=True)
    except TypeError:
        t, val = tref*np.nan, tref*np.nan
    return val


def longitude(pst):
    p = icp.get.params()
    t, val = icp.read.norm(pst, 'GPS', 'lon_ang', interp=True)
    return val

def latitude(pst):
    p = icp.get.params()
    t, val = icp.read.norm(pst, 'GPS', 'lat_ang', interp=True)
    return val

def roll(pst):
    p = icp.get.params()
    if pst.split('/')[1] in ['JKB2t']:
        if pst.split('/')[0] in ['SRH1', 'DEV', 'DEV2', 'HIC', 'NDEVON']:
            t, val = icp.read.targ(pst, 'treg', 'TRJ_JKB0', 'ztim_llzrphsaaa', interp=True, 
column=-7)
    else:
        t, val = icp.read.norm(pst, 'AVN', 'roll_ang', interp=True)
    return val

def signal(pst, pik, scale=1/1000., calib=True, air_loss=True, gain=0, **kwargs):
    """Extract signal from a pik file and apply various corrections
    """
    y, val = icp.read.pik(pst, pik, **kwargs)

    h = icp.get.surface_range(pst)
    # Pad the end of piks with nans to equal regular data length
    padding = np.abs((0,np.size(h)-np.size(val)))
    if len(h) > len(val):
        val = np.pad(val.astype(float), padding, 'constant', constant_values=np.nan)
    elif len(h) < len(val):
        h = np.pad(h.astype(float), padding, 'constant', constant_values=np.nan)
    else:
        pass
    print(scale)
    val = val*scale

    if calib is True:
        calval = 0.
        val = val +calval

    if air_loss is True:
        L = 10*np.log10( sr.utils.geo_loss(2*h) )
        gain = gain-L

    return val + gain


def ztim2frame(pst, year, day, ztim):
    """Convert ztim to the closest frame numbr in a PST
    """
    p = icp.get.params()
    z = icp.read.ztim(p['foc_path']+'/'+pst+'/ztim_DNhH')
    z_year  = z[1]
    z_day = z[3]
    z_ztim = z[5]
    w = [(z_year == year) * (z_day == day) * (z_ztim > ztim)]
    frame = np.nonzero(w)[1][0]
    return frame


#def surface_coefficients(pst, pik, wb=15e6, **kwargs):
#    """Surface coefficients (Reflectance and Scattering)
#    """
#    a = icp.read.rsr(pst, pik, **kwargs)
#    h = icp.get.surface_range(pst)[a['xo'].astype(int)]
#    Rsc, Rsn = rsr.invert.srf_coeff(Psc=a['pc'], Psn=a['pn'], h0=h, wb=15e6)
#    return {'Rsc':Rsc, 'Rsn':Rsn}


#def surface_properties(pst, pik, wf=60e6, **kwargs):
#    """Return surface permittivity and RMS height
#    """
#    a = icp.read.rsr(pst, pik, **kwargs)
#    h = icp.get.surface_range(pst)[a['xo'].astype(int)]
#    L = 10*np.log10( sr.utils.geo_loss(2*h) )
#    eps, sh = rsr.invert.spm(wf, a['pc']-L, a['pn']-L)
#    return {'sh':sh, 'eps':eps}

