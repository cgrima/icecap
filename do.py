"""
Various tools to process WAIS data
Author: Cyril Grima <cyril.grima@gmail.com>
"""

import fnmatch
import numpy as np
import icecap as icp
import os
import rsr.fit as fit
import rsr.utils as utils
import string
import subradar as sr
import time
import pandas as pd



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
       


def rsr(pst, pik, frame, **kwargs):
    """Apply RSR from a section of a transect
    """
    val = icp.get.signal(pst, pik, **kwargs)
    amp = 10**(val[frame[0]:frame[1]]/20.)
    return fit.lmfit(amp)



@loop
@timing
def rsr_inline(pst, pik, save=True, product='MagHiResInco1',**kwargs):
    """Process RSR along a profile
    """
    p = icp.get.params()

    val = icp.get.signal(pst, pik, product=product, **kwargs)
    amp = 10**(val/20)
    b = utils.inline_estim(amp, **kwargs)

    #Reorder data and remove geometric losses in air
    val2 = icp.get.signal(pst, pik, product=product, air_loss=False, **kwargs)
    diff = (val2-val)[b['xo']]

    data = pd.DataFrame({'1':b['xa'].values, 
                         '2':b['xo'].values,
                         '3':b['xb'].values, 
                         '4':b['pt'].values + diff,
                         '5':b['pc'].values + diff,
                         '6':b['pn'].values + diff,
                         '7':b['mu'].values, 
                         '8':b['crl'].values, 
                         '9':b['chisqr'].values,
                         '91':b['flag'].values})

    if save is True:
        folder = p['rsr_path'] + '/' + pst
        if not os.path.exists(folder):
            os.makedirs(folder)
        data.to_csv(folder + '/' + product+'.'+pik, sep='\t', float_format='%.7f',
                    na_rep='nan', header=False, index=False)
        print('CREATED: ' + folder + '/' + product+'.'+pik)



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
