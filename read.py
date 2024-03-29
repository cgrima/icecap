"""
Various tools to read WAIS data
Author: Cyril Grima <cyril.grima@gmail.com>
"""

import icecap as icp
import numpy as np
import os
import pandas as pd
#import string
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


def continuous_vec(x):
    """Modify a vector to make it continuous if needed. Designed to be
    applied on time stamps crossing midnight
    """
    i = np.argmin(x)
    if i != 0:
        x[i:] = x[i:] + x[i-1]
    return x


def norm(pst, instrument, stream, interp=False, **kwargs):
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
        data = np.interp(continuous_vec(foc_time), continuous_vec(time), data)
        time = foc_time

    return time, data


def tpro(pst, typ, fil, interp=True, **kwargs):
    """Read tpro or treg file
    !!!!!!!!!!!!!!!!!!
    !!! DEPRECATED !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Use read.targ instead !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    """
    p = icp.get.params()
    data_file = p['tpro_path'] + '/' + pst + '/' + typ + '/' + fil + '.bin'
    if icp.read.isfile(data_file) is False: return

    os.system('zvert ' + data_file + ' > tmp.tab')
    data = np.genfromtxt('tmp.tab')[:,-1]
    time = ztim('tmp.tab')['htim'].values
    os.system('rm tmp.tab')

    if interp is True: #interpolate to 1-m sampling
        foc_file = p['foc_path'] + '/' + pst + '/ztim_DNhH'

        if icp.read.isfile(foc_file) is False: return

        foc_time = ztim(foc_file)['htim'].values
        data = np.interp(continuous_vec(foc_time), continuous_vec(time), data)
        #data = np.interp(foc_time, time, data)
        time = foc_time

    return time, data


def targ(pst, folder, typ, fil, interp=True, column=-1, **kwargs):
    """Read a binary file from targ
    ARGUMENT examples
        pst : 'DEV/JKB2t/X101a'
        folder : 'treg'
        typ : 'TRJ_JKB0'
        fil : 'ztim_llzrphaaas'
    """
    p = icp.get.params()
    data_file = p[folder+'_path'] + '/' + pst + '/' + typ + '/' + fil + '.bin'
    if icp.read.isfile(data_file) is False: return

    os.system('zvert ' + data_file + ' > tmp.tab')
    data = np.genfromtxt('tmp.tab')[:,column]
    time = ztim('tmp.tab')['htim'].values
    os.system('rm tmp.tab')

    if interp is True: #interpolate to 1-m sampling
        foc_file = p['foc_path'] + '/' + pst + '/ztim_DNhH'

        if icp.read.isfile(foc_file) is False: return

        foc_time = ztim(foc_file)['htim'].values
        data = np.interp(continuous_vec(foc_time), continuous_vec(time), data)
        #data = np.interp(foc_time, time, data)
        time = foc_time

    return time, data


def ztim(fil):
    """Read time in a ztim-format file
    ARGUMENTS
        fil : string (full path + file name to read)
    """
    if icp.read.isfile(fil) is False: return
    out = pd.read_csv(fil, sep='\(|\)| |,', header=None)
    a = np.array(out)
    htim = (a[:,3] - a[:,3])*24 + a[:,5]*24/(86400*1e4)
    out['htim'] = htim
    return out


def pik(pst, pik, process=None, product='MagHiResInco1', **kwargs):
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
    fil = p['pik_path'].replace( p['process'], '') + process + '/' + pst + '/' + product + '.' + pik
    if icp.read.isfile(fil) is False: return
    os.system('grep P ' + fil + ' > ' + p['code_path'] + '/.read_pik.tmp')
    out = np.genfromtxt('.read_pik.tmp', delimiter='\t')
    os.system('rm ' + p['code_path'] + '/.read_pik.tmp')
    return out[:, 2], out[:, 3] # Y coordinate, value


def rsr(pst, pik, process=None, product='MagHiResInco1', **kwargs):
    """Read an rsr file
    """
    p = icp.get.params()
    if process is None:
        process = p['process']
    fil = p['rsr_path'].replace( p['process'], '') + process + '/' + pst + '/' + product + '.' + pik
    if icp.read.isfile(fil) is False: return

    a = np.genfromtxt(fil, delimiter='\t')
    out = {'xa':a[:,0], 'xo':a[:,1], 'xb':a[:,2], 'pt':a[:,3],'pc':a[:,4], 'pn':a[:,5], 'mu':a[:,6], 'crl':a[:,7], 'chisqr':a[:,8] }
    return out


def surface_coefficients(pst, pik, process=None, product='MagHiResInco1', **kwargs):
    """Read a surface_coefficients file
    """
    p = icp.get.params()
    if process is None:
        process = p['process']
    fil = p['rsr_path'].replace( p['process'], '') + process + '/' + pst + '/' + product + '.' + pik + '.surface_coefficients'
    if icp.read.isfile(fil) is False: return

    a = np.genfromtxt(fil, delimiter='\t')
    out = {'Rsc':a[:,0], 'Rsn':a[:,1], }
    return out


def surface_properties(pst, pik, process=None, product='MagHiResInco1', **kwargs):
    """Read a surface_coefficients file
    """
    p = icp.get.params()
    if process is None:
        process = p['process']
    fil = p['rsr_path'].replace( p['process'], '') + process + '/' + pst + '/' + product + '.' + pik + '.surface_properties'
    if icp.read.isfile(fil) is False: return
    a = np.genfromtxt(fil, delimiter='\t')
    out = {'sh':a[:,0], 'eps':a[:,1], 'flag':a[:,2]}
    return out


def bed_coefficients(pst, pik, process=None, product='MagHiResInco1', **kwargs):
    """Read a bed_coefficients file
    """
    p = icp.get.params()
    if process is None:
        process = p['process']
    fil = p['rsr_path'].replace( p['process'], '') + process + '/' + pst + '/' + product + '.' + pik + '.bed_coefficients'
    if icp.read.isfile(fil) is False: return

    a = np.genfromtxt(fil, delimiter='\t')
    out = {'Rbc':a[:,0], 'Rbn':a[:,1], }
    return out
