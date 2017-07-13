"""
Various tools to process WAIS data
Author: Cyril Grima <cyril.grima@gmail.com>
"""

import fnmatch
import numpy as np
import icecap as icp
import os
import rsr.fit as fit
import string
import time


def timing(func):
   """Outputs the time a function takes to execute.
   """
   def func_wrapper(*args, **kwargs):
       t1 = time.time()
       func(*args, **kwargs)
       t2 = time.time()
       print("- Processed in %.1f s.\n" % (t2-t1))
   return func_wrapper


def loop(func):
    """Decorator for processing over multiple PSTs
    """
    def func_wrapper(*args, **kwargs):
        pst_list = icp.get.pst(args[0])

        if 'from_process' in locals():
            process = from_process
        elif 'process' not in locals():
            process = None
            for pst_i in pst_list:
                product_list, pik_list = icp.get.pik(pst_i, process=process, **kwargs)
                for i, pik_i in enumerate(pik_list):
                    if fnmatch.filter([pik_i], args[1]):
                        func(pst_i, pik_i, **kwargs)

    return func_wrapper


def signal_calibration(val, scale=1/1000., db=True, power=True):
    """Calibrates the signal from pik files (Energy in dB) with predifined values.
    ARGUMENTS
        val : float or list of floats (Energy values)
        dB : bool (if True, in dB, if False, linear values)
        power : bool (if True, power, if False, amplitude values)
    OUTPUT
       Calibrated Energy in dB (Does not include geometric spreading) 
    """
    out = val*scale
    if power is True:
        if db is True:
           out = out
        elif db is False:
           out = 10**(out)
    elif power is False:
        if db is True:
            out = out/20.
        elif db is False:
            out = 10**(out/20.)

    return out


#@timing
def rsr(pst, pik, lim, calib=True, gain=0, fit_model='hk', bins='knuth', **kwargs):
    """Apply RSR from a section of a transect
    ARGUMENTS
        calval : bool (Tell to use the calibration value stored in the hierarchy)
        gain : float (Energy in dB to be added to the data)
    """
    y, val = icp.read.pik(pst, pik, **kwargs)
    if calib is True:
        val = icp.do.signal_calibration(val)

    amp = 10**(val/20.+gain/20.)
    a = fit.lmfit(amp[lim[0]:lim[1]], fit_model=fit_model, bins=bins)
    return a


@loop
@timing
def rsr_inline(pst, pik, process=None, product='MagHiResInco1', save=True, winsize=1000., sampling=250., verbose=True, **rsr_kwargs):
    """Process RSR along a profile
    """
    y, val = icp.read.pik(pst, pik, process=process, product=product)
    amp = icp.do.signal_calibration(val, db=False, power=False)
    return amp


@loop
@timing
def topik1m(pst, pik, from_process='pik1.RADnh3', from_product='MagLoResInco1', to_product='MagHiResInco1'):
    """ Inteprolate any pik file into 1m sampling
    """
    p = icp.get.params()

    source = string.replace(p['pik_path'], p['process'], '') + from_process + '/' + pst + '/' + from_product + '.' + pik
    bxds = p['cmp_path'] + '/' + pst + '/'+ to_product

    if icp.read.isfile(source) is False: return
    if icp.read.isfile(bxds) is False: return

    if not os.path.exists(p['pik_path'] + '/' + pst):
        os.makedirs(p['pik_path'] + '/' + pst)

    target = p['pik_path'] + '/' + pst + '/'+ to_product + '.' + pik
    LU = target + '_LU'
    P = target + '_P'
    os.system('mkdir -p ' + p['pik_path'] + '/' + pst)
    os.system('pik4Hzto1m ' + pst + ' < ' + source + ' > ' + LU)
    os.system('pk3 3200 0 3200 ' + bxds + ' < ' + LU + ' > ' + P)
    os.system('cat ' + LU + ' ' + P + ' > ' + target)
    os.system('rm ' + LU + ' ' + P)
    print('CREATED: ' + target)
