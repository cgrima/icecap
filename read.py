"""
Various tools to read WAIS data
Author: Cyril Grima <cyril.grima@gmail.com>
"""

import icecap as icp
import numpy as np
import os
import pandas as pd
import string
import sys


def isfile(fil, verbose=True):
    """exit the process if file does not exist
    """
    out = os.path.isfile(fil)
    if out is False:
        message = fil + ' does not exist'
        if verbose is True:
            print('IGNORED: ' + message)
    return out


def norm(pst, instrument, stream, interp=False):
    """Read data streams from norm
    ARGUMENTS
        pst : string (e.g. 'MIS/JKB2e/Y35a')
        instrument : string (e.g. 'LAS')
        stream : string (e.g. 'las_rng')
        1m : bool (True is for interpolation to 1m sampling)
    """
    p = icp.get.params()

    data_file = p['norm_path'] + '/' + pst + '/' + instrument + '_' + \
                pst.split('/')[1][0:3] + 'a/' + stream
    time_file = p['norm_path'] + '/' + pst + '/' + instrument + '_' + \
               pst.split('/')[1][0:3] + 'a/syn_ztim'

    if icp.read.isfile(data_file) is False: return
    if icp.read.isfile(time_file) is False: return

    data = np.genfromtxt(data_file)
    time = ztim(time_file)['htim'].values

    if interp is True: #interpolate to 1-m sampling
        foc_file = p['foc_path'] + '/' + pst + '/ztim_DNhH'

        if icp.read.isfile(foc_file) is False: return
   
        foc_time = ztim(foc_file)['htim'].values
        data = np.interp(foc_time, time, data)
        time = foc_time

    return time, data


def ztim(fil):
    """Read time in a ztim-format file
    ARGUMENTS
        fil : string (full path + file name to read)
    """
    if icp.read.isfile(fil) is False: return
    out = pd.read_csv(fil, sep='\(|\)| |,', header=None, engine='python')
    a = np.array(out)
    htim = (a[:,3] - a[:,3])*24 + a[:,5]*24/(86400*1e4)
    out['htim'] = htim
    return out


def pik(pst, pik, process=None, product='MagHiResInco1'):
    """Read pick files
    ARGUMENTS
        pst : string (e.g. 'MIS/JKB2e/Y35a')
        pik : string (e.g. 'srf_elg')
        process : string (e.g. 'pik1.1m.RADnh3')
        product : string (e.g. 'MagHiResInco1')
    OUTPUT
        list : list (Y coordinate, value) 
    """
    p = icp.get.params()
    if process is None:
        process = p['process']
    fil = string.replace(p['pik_path'], p['process'], '') + process + '/' + pst + '/' + product + '.' + pik
    if icp.read.isfile(fil) is False: return
    os.system('grep P ' + fil + ' > ' + p['code_path'] + '/.read_pik.tmp')
    out = np.genfromtxt('.read_pik.tmp', delimiter='\t')
    os.system('rm ' + p['code_path'] + '/.read_pik.tmp')
    return out[:, 2], out[:, 3] # Y coordinate, value
