"""
Various tools to manipulate WAIS data
Author: Cyril Grima <cyril.grima@gmail.com>
"""

import numpy as np
import pandas as pd
import os
import glob
import string
import rsr.utils
import raw
import scipy.constants as ct
import matplotlib.pyplot as plt
from snow import properties, z_profile
from params import *


def bathymetry(pst, ext1, ext2, stat='hk', inv='spm', save=True):
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
    stat : string
        statistical method used
    inv : string
        backscattering model used for physical properties inversion
    save : bool
        wether to save the results in a txt file
    """
    rsr1 = raw.read_rsr(pst, ext1, stat=stat, inv=inv)
#    rsr2 = raw.read_rsr(pst, ext2, stat=stat, inv=inv)
    z1, pik1 = raw.read_pik(pst, ext1)
    z2, pik2 = raw.read_pik(pst, ext2, process='MagHiResInco1')

    lat = rsr1.lat#.values.tolist()
    lon = rsr1.lon#.values.tolist()
    eps = rsr1.eps#.values.tolist()
    crl = rsr1.crl#.values.tolist()
    sh = rsr1.eps#.values.tolist()
    roll = rsr1.roll#.values.tolist()
    flag = rsr1.flag#.values.tolist()
    w = rsr1.xa.values.tolist()
    delay = abs(z2[w]-z1[w])*20e-9
    delay[np.isnan(delay)] = 0.
    uncertainty = eps*0.
    depth = eps*0.
    dns_at_depth = eps*0.
    
    # Depth conversion
    for i, val in enumerate(eps):
        if np.isnan(val) is False:
            z = np.linspace(0, 5000, 5000)
            dns0 = properties.perm_kovacs95(val)
            dns_a = z_profile.density_sorgelaw(dns0, z, 10./1.9)
            dns_b = z_profile.density_sorgelaw(dns0, z, 68./1.9)
            eps_a = properties.perm_kovacs95(dns_a, density=True)
            eps_b = properties.perm_kovacs95(dns_b, density=True)
            dt_a = np.cumsum(2*np.sqrt(eps_a)/ct.c)
            dt_b = np.cumsum(2*np.sqrt(eps_b)/ct.c)
            w_a = np.where(dt_a >= delay[i])[0][0]
            w_b = np.where(dt_b >= delay[i])[0][0]
            depth_a = z[w_a]
            depth_b = z[w_b]
            uncertainty[i] = abs(depth_a-depth_b)/2. +  ct.c/(2*bdw*np.sqrt(eps[i]))/2.
            depth[i] = np.min([depth_a, depth_b])+ abs(depth_a-depth_b)/2.
            dns_at_depth[i] = np.mean([dns_a[w_a], dns_b[w_b]])
    
    # Results
    out = {'lat':lat, 'lon':lon, 'eps':eps, 'sh':sh, 'crl':crl, 'delay':delay, 
           'uncertainty':uncertainty, 'depth':depth, 'dns_at_depth':dns_at_depth,
           'roll':roll, 'flag':flag}
    out = pd.DataFrame(out)

    if save is True:
        save_fil = rsr_path + '/' + string.join(pst.split('/'), '_') + '.' + \
                   ext1 + '-' + ext2
        out.to_csv(save_fil + '.bthm.txt', sep='\t', index=False, float_format='%.7f')

    return out


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


def inline_rsr(pst, ext, stat='hk', inv='spm', save=True, winsize=1000.,
               sampling=250., verbose=True, process=process, other_gain=logain, **kwargs):
    """launch sliding RSR along a track

    Arguments
    ---------
    pst : string
        pst name (e.g. 'MIS/JKB2e/Y35a')
    ext : string
        pik file extension (e.g. 'elg_brn')
    
    Keywords
    --------
    save : bool
        wether or not to save the results
    """
    geo = raw.read_geo(pst)
    y, val = raw.read_pik(pst, ext, process=process)
    rng = geo.rng[np.arange(val.size)]
    amp = 10**(calibration(val, rng=rng, other=other_gain)/20.)
    b = rsr.utils.inline_estim(amp, frq=frq, stat=stat, inv=inv, winsize=winsize, sampling=sampling, verbose=verbose, **kwargs)
    xo = np.round(np.array(b.xo)) # positions of the computed statistics
    b['lat'] = np.array(geo.ix[xo, 'lat'])
    b['lon'] = np.array(geo.ix[xo, 'lon'])
    b['roll'] = np.array(geo.ix[xo, 'roll'])
    b['rng'] = np.array(geo.ix[xo, 'rng'])

    if save is True:
        save_fil = string.replace(os.getcwd(), 'code', 'targ') + '/' + \
                   string.join(pst.split('/'), '_') + '.' + ext + '.' + stat + '.' + inv
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

