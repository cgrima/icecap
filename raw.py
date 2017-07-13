"""

!!!!!!!!!!!!!!!!
!! DEPRECATED !!
!!!!!!!!!!!!!!!!

Various tools to read WAIS data
Author: Cyril Grima <cyril.grima@gmail.com>
"""

import numpy as np
import pandas as pd
import os
import glob
import string
from scipy.interpolate import UnivariateSpline
from params import *


def grep_params(word, season = season):
    """grep a line in 'params' file for a given season
    
    Arguments
    ---------
    word : string
        grep line that include 'word'

    Keywords
    --------
    season : string
        season name (e.g. 'ICP4') 
    """
    a = code_path.split('/')
    fil = string.join(a[0: -5], '/') + '/pcor' + '/' + season + '/params'
    out = os.popen('grep ' + word + ' ' + fil).read()
    out = string.replace(out, '\n', '')
    out = string.replace(out, '#', '')
    out = string.replace(out, '(', '')
    out = string.replace(out, ')', '')
    out = string.replace(out, '{', '')
    out = string.replace(out, '}', '')
    out = string.replace(out, ',', '')
    return out.split()


def list_pst(folder, pst, ext='*'):
    """list PSTs present in various folders

    Arguments
    ---------
    folder : string
        folder name not including $WAIS (e.g. 'orig/xtra/ICP4/PIK/pik1.1m.RADnh3')
    pst : string
        pst name (e.g. 'MIS/JKB2e/Y35a')

    Keywords
    --------
    ext : string
        file extension if needed
    """
    if folder.find('/RSR/') is not -1:
        template = string.join([wais_path, folder, ''], '/') + \
                   string.join(pst.split('/'), '_') + '.' + ext + '*'
    else:
        template = string.join([wais_path, folder, pst, ''], '/') + '*' + ext
    print(template)    
    return glob.glob(template)


def read_ztim(fil):
    """Read time in a ztim-format file with the following structure:
    
    Arguments
    ---------
    fil : string
        full path + file name to read
    """
    out = pd.read_csv(fil, sep='\(|\)| |,', header=None)
    a = np.array(out)
    htim = (a[:,3] - a[:,3])*24 + a[:,5]*24/(86400*1e4)
    out['htim'] = htim
    return out


def read_geo(pst):
    """read some geographic parameters from nrm and interpolate them
    with the 1-m data ztim

    Arguments
    ---------
    pst : string
        pst name (e.g. 'MIS/JKB2e/Y35a')
    """

    AC = pst.split('/')[1][0:3]

    foc_time = read_ztim(foc_path + '/' + pst + '/ztim_DNhH').htim
    avn_time = read_ztim(norm_path + '/' + pst + '/AVN_'+AC+'a/syn_ztim').htim
    las_time = read_ztim(norm_path + '/' + pst + '/LAS_'+AC+'a/syn_ztim').htim
    
    lat_nrm = np.genfromtxt(norm_path + '/' + pst + '/AVN_'+AC+'a/lat_ang')
    lon_nrm = np.genfromtxt(norm_path + '/' + pst + '/AVN_'+AC+'a/lon_ang')
    roll_nrm = np.genfromtxt(norm_path + '/' + pst + '/AVN_'+AC+'a/roll_ang')
    rng_nrm = np.genfromtxt(norm_path + '/' + pst + '/LAS_'+AC+'a/las_rng')

    lat = np.interp(foc_time, avn_time, lat_nrm)
    lon = np.interp(foc_time, avn_time, lon_nrm)
    roll = np.interp(foc_time, avn_time, roll_nrm)
    rng = np.interp(foc_time, las_time, rng_nrm)
    #lat = UnivariateSpline(avn_time, lat_nrm)(foc_time)
    #lon = UnivariateSpline(avn_time, lon_nrm)(foc_time)
    #roll = UnivariateSpline(avn_time, roll_nrm)(foc_time)
    #rng = UnivariateSpline(las_time, rng_nrm)(foc_time)

    return pd.DataFrame({'lat':lat, 'lon':lon, 'roll':roll, 'rng':rng})


def read_pik(pst, ext, process='MagHiResInco1'):
    """Read pick files

    Arguments
    ---------
    pst : string
        pst name (e.g. 'MIS/JKB2e/Y35a')
    ext : string
        pik file extension (e.g. 'elg_brn')

    Keywords
    --------
    process : string
        process name

    Output
    ------
    list, list
        Y coordinate, value 
    """
    fil = pik_path + '/' + pst + '/' + process + '.' + ext
    os.system('grep P ' + fil + ' > ' + code_path + '/.read_pik.tmp')
    b = np.genfromtxt('.read_pik.tmp', delimiter='\t')
    return b[:, 2], b[:, 3] # Y coordinate, value


def read_rsr(pst, ext, fit_model='hk', inv='spm'):
    """Read RSR files

    Arguments
    ---------
    pst : string
        pst name (e.g. 'MIS/JKB2e/Y35a')
    ext : string
        pik file extension (e.g. 'elg_brn')

    Keywords
    --------
    fit_model : string
        statistical method used
    inv : string
        backscattering model used for physical properties inversion 
    """
    fil = rsr_path+'/'+string.join(pst.split('/'), '_') + '.' + ext + '.' + fit_model + '.' + inv + '.txt'
    out = pd.read_table(fil)
    return out


def read_bthm(pst, ext1, ext2):
    """Read BTHM files

    Arguments
    ---------
    pst : string
          pst name (e.g. 'MIS/JKB2e/Y35a')
    ext1 : string
        pik file extension (e.g. 'elg_brn') for the 1st arrival
    ext2 : string
        pik file extension (e.g. 'elg_brn') for the 2nd arrival
    """
    fil = glob.glob(rsr_path+'/'+string.join(pst.split('/'), '_') + '.' + ext1 + '-' + \
          ext2 + '.bthm.txt')[0]
    out = pd.read_table(fil)
    return out
