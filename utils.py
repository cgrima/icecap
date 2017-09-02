"""

!!!!!!!!!!!!!!
! DEPRECATED !
!!!!!!!!!!!!!!

Various tools to manipulate WAIS data
Author: Cyril Grima <cyril.grima@gmail.com>
"""

import numpy as np
import pandas as pd
import os
import glob
import string
import rsr.utils
from icecap import raw
import scipy.constants as ct
import matplotlib.pyplot as plt
from snow import properties, z_profile
from params import *


def bathymetry(pst, ext1, ext2, fit_model='hk', inv='spm', save=True):
    """Compute bathymetry

    Arguments
    ---------
    pst : string
        pst pattern (example: 'MIS_JKB2e*')
    ext1 : string
        file extension for the top interface (e.g. 'srf_cyg.hk.spm')
    ext2 : string
        file extension for the bottom interface (e.g. 'srf_cyg.hk.spm')

    Keywords
    --------
    fit_model : string
        statistical method used
    inv : string
        backscattering model used for physical properties inversion
    save : bool
        wether to save the results in a txt file
    """
    rsr1 = raw.read_rsr(pst, ext1, fit_model=fit_model, inv=inv)
    z1, pik1 = raw.read_pik(pst, ext1)
    z2, pik2 = raw.read_pik(pst, ext2, process='MagHiResInco1')

    # Parameters
    lat = rsr1.lat
    lon = rsr1.lon
    eps = rsr1.eps
    crl = rsr1.crl
    sh = rsr1.eps
    roll = rsr1.roll
    flag = rsr1.flag
    w = rsr1.xa.values.tolist()
    delay = abs(z2[w]-z1[w])*20e-9
    delay[np.isnan(delay)] = 0.

    # Results    
    z = [dns_depth(eps[i], delay[i]) for i, val in enumerate(eps)]
    depth = [i[0] for i in z]
    uncertainty = [i[1] for i in z]
    dns_at_depth = [i[2] for i in z]

    out = {'lat':lat, 'lon':lon, 'eps':eps, 'sh':sh, 'crl':crl, 'delay':delay, 
           'uncertainty':uncertainty, 'depth':depth, 'dns_at_depth':dns_at_depth,
           'roll':roll, 'flag':flag}
    out = pd.DataFrame(out)

    if save is True:
        save_fil = rsr_path + '/' + string.join(pst.split('/'), '_') + '.' + \
                   ext1 + '-' + ext2
        out.to_csv(save_fil + '.bthm.txt', sep='\t', index=False, float_format='%.7f')

    return out


def dns_depth(eps, delay, z=np.linspace(0, 5000, 5000), zp=[19., 129.]):
    """Give an estimation for the depth of a subsurface echo and based on the surface density.
    Assumes a Sorge's law depth/density profile

    Arguments
    ---------
    eps : float
        permittivity at the surface (z=0)
    delay : Float
        time delay between the surface and subssurface echo [sec.]

    Keywords
    --------
    z : array of float
        Depth profile to consider [m]
    zp : [Float, Float]
        Depth of the firn/ice transition to consider for the Sorge's Law

    Output
    ------
    depth : Estimated depth for the subsurface echo [m]
    uncertainty: +/- uncertainty on depth
    dns_at_depth: Estimated density at depth [kg.m^{-3}] 
    """
    if np.isfinite(eps) == False: eps = 0.
    if np.isfinite(delay) == False: delay = 0.
    # Density conversions
    dns0 = properties.perm_kovacs95(eps)
    dns_a = z_profile.density_sorgelaw(dns0, z, zp[0])
    dns_b = z_profile.density_sorgelaw(dns0, z, zp[1])
    # Permittivity conversions
    eps_a = properties.perm_kovacs95(dns_a, density=True)
    eps_b = properties.perm_kovacs95(dns_b, density=True)
    # Time profiles
    dt_a = np.cumsum(2*np.sqrt(eps_a)/ct.c)
    dt_b = np.cumsum(2*np.sqrt(eps_b)/ct.c)
    # Depth estimations
    w_a = np.where(dt_a >= delay)[0][0]
    w_b = np.where(dt_b >= delay)[0][0]
    depth_a = z[w_a]
    depth_b = z[w_b]
    
    uncertainty = abs(depth_a-depth_b)/2. +  ct.c/(2*bdw*np.sqrt(eps))/2.
    depth = np.min([depth_a, depth_b])+ abs(depth_a-depth_b)/2.
    dns_at_depth = np.mean([dns_a[w_a], dns_b[w_b]])
    
    return depth, uncertainty, dns_at_depth


def calibration(val, scale=scale, wl=ct.c/frq, pt=pt, antenna=antgain*2, rng = False,
                abs_calib = abs_calib, other=logain):
    """Signal calibrated from instrumental and geophysic gains (power in dB)
    If rng is given, will correct for the 2-way specular geomtric losses

    Arguments
    ---------
    val : float or list of float
        raw echo value(s)

    Keywords
    --------
    scale : float
        scale factor
    wl : float
        wavelength [m]
    pt : float
        transmitted power [dB]
    antenna : float
        antenna gain [dB]
    rng : float
        Range to the target [m]
    abs_calib : float
        Absolute calibration value usually obtained over a known terrain [dB]
    other : float
        other gains [dB]
    """
    geometric_loss = 0 if rng is False else 20*np.log10(8*np.pi*rng)
    out = val*scale - 20*np.log10(wl) - pt - antenna + geometric_loss + abs_calib + other
    return out


def dns2eps(dns):
    """Convert dry-snow density into permittivity (epsilon)
    Based on Kovacs et al. [1995]
    
    Arguments
    ---------
    dns : float
        density [kg.m^{-3}]
    """
    return (1+845e-6*dns)**2


def eps2dns(eps):
    """Convert permittivity (epsilon) into dry-snow density [kg.m^{-3}] 
    Based on Kovacs et al. [1995]
    
    Arguments
    ---------
    eps : float
        density [kg.m^{-3}]
    """
    return (np.sqrt(eps)-1)/845e-6


def inline_rsr(pst, ext, fit_model='hk', inv='spm', save=True, winsize=1000.,
               sampling=250., verbose=True, process='MagHiResInco1', other_gain=logain,
               rng='surface', **kwargs):
    """launch sliding RSR along a track

    Arguments
    ---------
    pst : string
        pst name (e.g. 'MIS/JKB2e/Y35a')
    ext : string
        pik file extension (e.g. 'elg_brn') for the echo
    
    Keywords
    --------
    save : bool
        wether or not to save the results
    """
    #--------------------------------------------------------------------------
    # Data Retrieval
    #--------------------------------------------------------------------------
    geo = raw.read_geo(pst)
    y, val = raw.read_pik(pst, ext, process=process)
    if rng is 'surface':    
        rng = geo.rng[np.arange(val.size)]
    if rng is 'bed':
        bthm = raw.read_bthm(pst, '*', ext)
        air_rng = geo.rng[np.arange(val.size)]
        x, xp, fp = np.arange(val.size), np.arange(bthm.shape[0]), np.array(bthm.depth)
        sub_rng = np.interp(x, xp, fp)
        rng = air_rng + sub_rng
    amp = 10**(calibration(val, rng=rng, other=other_gain)/20.)
    #--------------------------------------------------------------------------
    # Data Retrieval
    #--------------------------------------------------------------------------
    b = rsr.utils.inline_estim(amp, frq=frq, fit_model=fit_model, inv=inv, winsize=winsize, sampling=sampling, verbose=verbose, **kwargs)
    xo = np.round(np.array(b.xo)) # positions of the computed statistics
    b['lat'] = np.array(geo.ix[xo, 'lat'])
    b['lon'] = np.array(geo.ix[xo, 'lon'])
    b['roll'] = np.array(geo.ix[xo, 'roll'])
    b['rng'] = np.array(geo.ix[xo, 'rng'])

    if save is True:
        save_fil = string.replace(os.getcwd(), 'code', 'targ') + '/' + \
                   string.join(pst.split('/'), '_') + '.' + ext + '.' + fit_model + '.' + inv
	title = pst + ' ' + ext
        b.to_csv(save_fil + '.txt', sep='\t', index=False, float_format='%.7f')
        rsr.utils.plot_inline(b, frq=frq, title=title)
        plt.savefig(save_fil + '.png', bbox_inches='tight')
    return b


def group(pst, ext, save=True, rem_bad=True):
    """Group  data of a given type in one txt file

    Arguments
    ---------
    pst : string
        pst pattern (example: 'MIS_JKB2e*')
    ext : string/disk/kea/WAIS/orig/xtra/ICP4/PIK/pik1.1m.RADnh3/MIS/JKB2e/Y4
        file extension (example: 'srf_cyg.hk.spm')

    Keywords
    --------
    save : bool
        save the data
    rem_bad : bool
        remove rows with bad values
    """
    files = glob.glob(rsr_path+'/'+string.replace(pst, '/', '_')+'.'+ext+'.txt')
    out = pd.read_table(files[0])
    for fil in files:
        if fil != files[0]:
            a = pd.read_table(fil)
            out = pd.concat([out, a])
    if rem_bad is True:
        out = out.replace([np.inf, -np.inf], np.nan)
        out = out.dropna()
	out = out[np.abs(out.roll) < 2]
	out = out[out.flag == 1]
    if save is True:
        name = string.replace(string.replace(pst, '/', '_'), '*', '').rstrip('_')
        filename = rsr_path+'/'+name+'.'+ext+'.txt'
        out.to_csv(filename, sep='\t', index=False, float_format='%.7f')
    return out


def topik1m(pst, ext, process0='pik1.RADnh3', process='MagLoResInco1', process_target='MagHiResInco1'):
    """Group  data of a given type in one txt file

    Arguments
    ---------
    pst : string
        pst pattern (example: 'MIS_JKB2e*')
    ext : string/disk/kea/WAIS/orig/xtra/ICP4/PIK/pik1.1m.RADnh3/MIS/JKB2e/Y4
        file extension (example: 'srf_cyg.hk.spm')
    process0 ; string
        process0 of the pik file to convert
    process : string
        process of the pik file to convert
    
    Example
    -------
    topik1(pst, 'srf_elg', process0='pyk1.RADnh3', process='MagLoResInco1')
    """
    source = os.path.split(pik_path)[0] + '/' + process0 + '/' + pst + '/' + process + '.' + ext
    bxds = cmp_path + '/' + pst + '/'+process_target
    target = pik_path + '/' + pst + '/'+process_target+'.' + ext
    LU = target + '_LU'
    P = target + '_P'

    os.system('mkdir -p ' + pik_path + '/' + pst)
    os.system('pik4Hzto1m ' + pst + ' < ' + source + ' > ' + LU)
    os.system('pk3 3200 0 3200 ' + bxds + ' < ' + LU + ' > ' + P)
    os.system('cat ' + LU + ' ' + P + ' > ' + target)
    os.system('rm ' + LU + ' ' + P)
