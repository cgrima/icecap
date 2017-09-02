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

        if 'from_process' in kwargs:
            process = kwargs['from_process']
        elif 'process' not in kwargs:
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
def topik1m(pst, pik, from_process='pik1', from_product='MagLoResInco1', to_product='MagHiResInco1'):
    """ Inteprolate any pik file into 1m sampling
        NOTE: at this point, from_process is mandatory for @loop to work
    """
    p = icp.get.params()

    source = string.replace(p['pik_path'], p['process'], '') + from_process + '/' + pst + '/' + from_product + '.' + pik
    bxds = p['cmp_path'] + '/' + pst + '/'+ to_product

    test = icp.read.isfile(source) * icp.read.isfile(bxds)
    if test is 0: return

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



@loop
@timing
def gather(pst, pik, fil=None, product='MagHiResInco1', **kwargs):
    """Gather data in a text file
    """
    p = icp.get.params()
    a = pd.DataFrame()

    print('Gathering data from '+pst)

    r = icp.read.rsr(pst, pik, **kwargs)
    if os.path.isfile(p['rsr_path'] + '/' + pst + '/' + product + '.' + pik ):
      xo = [np.int(i) for i in r['xo']]
      #a['pst'] = np.array([pst for i in xo])
      a['xo'] = xo
      a['crl'] = r['crl']
      a['longitude'] = icp.get.longitude(pst)[xo]
      a['latitude'] = icp.get.latitude(pst)[xo]
      a['roll'] = icp.get.roll(pst)[xo]
      a['surface_range'] = icp.get.surface_range(pst)[xo]
    else:
      return
    
    if os.path.isfile(p['rsr_path'] + '/' + pst + '/' + product + '.' + pik + '.surface_coefficients'):
        b = icp.read.surface_coefficients(pst, pik, **kwargs)
        a['Rsc'] = b['Rsc']
        a['Rsn'] = b['Rsn']
    else:
        a['Rsc'] = xo*np.nan
        a['Rsn'] = xo*np.nan

    if os.path.isfile(p['rsr_path'] + '/' + pst + '/' + product + '.' + pik + '.surface_properties'):
        b = icp.read.surface_properties(pst, pik, **kwargs)
        a['sh'] = b['sh']
        a['eps'] = b['eps']
        a['flag'] = r['flag']*b['flag']
    else:
        a['sh'] = xo*np.nan
        a['eps'] = xo*np.nan
        a['flag'] = xo*np.nan

    if fil is None:
        fil = p['season'] + '_gather.csv'

    if os.path.isfile(fil):
        c = pd.read_csv(fil)
        d = c.append(a, ignore_index=True)
        d.to_csv(fil, index=False, float_format='%.7f', na_rep='nan')
    else:
        a.to_csv(fil, index=False, float_format='%.7f', na_rep='nan')
  


@loop
@timing
def surface_coefficients(pst, pik, wb=15e6, product='MagHiResInco1', save=True, **kwargs):
    """Surface coefficients (Reflectance and Scattering)
    """
    p = icp.get.params()
    a = icp.read.rsr(pst, pik, product='MagHiResInco1', **kwargs)
    h = icp.get.surface_range(pst)[a['xo'].astype(int)]
    Rsc, Rsn = invert.srf_coeff(Psc=a['pc'], Psn=a['pn'], h0=h, wb=15e6)
 
    out = {'0_Rsc':Rsc, '1_Rsn':Rsn}

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

    eps, sh, flag = np.zeros(np.size(h)), np.zeros(np.size(h)), np.zeros(np.size(h))

    for i, val in enumerate(h):
        tmp = invert.spm(wf, a['pc'][i]-L[i], a['pn'][i]-L[i])
        eps[i] = tmp['eps']
        sh[i] = tmp['sh']
        flag[i] = np.int(3e8/wf*0.05 > tmp['sh'])

    out = {'0_sh':sh, '1_eps':eps, '2_flag':flag}

    if save is True:
        p = icp.get.params()
        target = p['rsr_path'] + '/' + pst + '/' + product + '.' + pik + '.' + inspect.stack()[0][3]
        icp.do.save(out, target)



@loop
@timing
def bed_coefficients(pst, bed_pik, wb=15e6, product='MagHiResInco1', save=True, **kwargs):
    p = icp.get.params()
    # get first available srf pik
    foo, pik_list = icp.get.pik('SMIS/MKB2l/Y51a')
    srf_pik_list = [i for i in bar if 'srf' in i]

    if srf_pik_list:
        srf_pik = srf_pik_list[0]
    else:
        print('IGNORED: No srf pik for '+pst)
        return

        s = icp.read.rsr(pst, srf_pik, product='MagHiResInco1', **kwargs)
    b = icp.read.rsr(pst, bed_pik, product='MagHiResInco1', **kwargs)

    h0 = icp.get.surface_range(pst)[a['xo'].astype(int)]
    #h1 = 
    #n1 = 
    #sh = 
    #Q1 = 
    
    Rbc, Rbn = invert.bed_coeff(Psc=s['pc'], Psn=s['pn'],
                                Pbc=b['pc'], Pbn=b['pn'],
                                n1=n1, sh=sh, h0=h0, h1=h1, Q1=Q1,
                                wb=15e6)

    out = {'0_Rbc':Rsc, '1_Rbn':Rsn}
