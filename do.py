"""
Various tools to process WAIS data
Author: Cyril Grima <cyril.grima@gmail.com>
"""

import fnmatch
import numpy as np
import icecap as icp
import inspect
import os
import rsr.fit as fit
import rsr.utils as utils
import rsr.invert as invert
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




def save(data, target):
    """Save a dictionnary of data in a text file
    """
    df = pd.DataFrame(data)
    df.to_csv(target, sep='\t', index=False, float_format='%.7f', header=False, na_rep='nan')
    print('CREATED: ' + target)
    
                       


def rsr(pst, pik, frame, **kwargs):
    """Apply RSR from a section of a transect
    """
    val = icp.get.signal(pst, pik, **kwargs)
    amp = 10**(val[frame[0]:frame[1]]/20.)
    return fit.lmfit(amp)



@loop
@timing
def rsr_inline(pst, pik, save=True, product='MagHiResInco1', **kwargs):
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


def gather(pst, pik, file=None, **kwargs):
    """Gather data in a text file
    """
    p = icp.get.params()
    a = pd.DataFrame()

    if os.path.isfile():
       pass 
    #if os.path.isfile(file) is False:
    #else 
  


@loop
@timing
def surface_coefficients(pst, pik, wb=15e6, product='MagHiResInco1', save=True, **kwargs):
    """Surface coefficients (Reflectance and Scattering)
    """
    p = icp.get.params()
    a = icp.read.rsr(pst, pik, product='MagHiResInco1', **kwargs)
    h = icp.get.surface_range(pst)[a['xo'].astype(int)]
    Rsc, Rsn = invert.srf_coeff(Psc=a['pc'], Psn=a['pn'], h0=h, wb=15e6)
    out = {'Rsc':Rsc, 'Rsn':Rsn}

    if save is True:
        p = icp.get.params()
        target = p['rsr_path'] + '/' + pst + '/' + product + '.' + pik + '.' + inspect.stack()[0][3]
        icp.do.save(out, target)



@loop
@timing
def surface_properties(pst, pik, wf=60e6, product='MagHiResInco1', save=True, **kwargs):
    """Return surface permittivity and RMS height
    """
    a = icp.read.rsr(pst, pik, product='MagHiResInco1', **kwargs)
    h = icp.get.surface_range(pst)[a['xo'].astype(int)]
    L = 10*np.log10( sr.utils.geo_loss(2*h) )

    eps, sh = np.zeros(np.size(h)), np.zeros(np.size(h))

    for i, val in enumerate(h):
        tmp = invert.spm(wf, a['pc'][i]-L[i], a['pn'][i]-L[i])
        eps[i] = tmp['eps']
        sh[i] = tmp['sh']

    out = {'sh':sh, 'eps':eps}

    if save is True:
        p = icp.get.params()
        target = p['rsr_path'] + '/' + pst + '/' + product + '.' + pik + '.' + inspect.stack()[0][3]
        icp.do.save(out, target)

