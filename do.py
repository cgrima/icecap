"""
Various tools to process WAIS data
Author: Cyril Grima <cyril.grima@gmail.com>
"""

import fnmatch
import numpy as np
import icecap as icp
import inspect
import os
import rsr.run as run
import rsr.fit as fit
import rsr.utils as utils
import rsr.invert as invert
#import string
import subradar as sr
import time
import pandas as pd
import multiprocessing



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
    """Decorator for processing sequentialy over multiple PSTs
    """
    def func_wrapper(*args, **kwargs):
        pst_list = icp.get.pst(args[0])

        if 'from_process' in kwargs:
            process = kwargs['from_process']
        elif 'process' not in kwargs:
            process = None

        _processes = []
        for pst_i in pst_list:
            product_list, pik_list = icp.get.pik(pst_i, process=process, **kwargs)
            for i, pik_i in enumerate(pik_list):
                if fnmatch.filter([pik_i], args[1]):
                    func(pst_i, pik_i, **kwargs)
#                    p = multiprocessing.Process(target=func, args=(pst_i, pik_i), kwargs=kwargs)
#                    p.start()
#                    _processes.append(p)
#        for p in _processes:
#            p.join()
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
    return run.processor(amp)
    #return fit.lmfit(amp)



@loop
@timing
def rsr_inline(pst, pik, save=True, product='MagHiResInco1', nbcores=4, **kwargs):
    """Process RSR along a profile
    """
    p = icp.get.params()

    val = icp.get.signal(pst, pik, product=product, air_loss=False, **kwargs)
    amp = 10**(val/20)
    b = run.along(amp, nbcores=nbcores, **kwargs)

    data = pd.DataFrame({'1':b['xa'].values,
                         '2':b['xo'].values,
                         '3':b['xb'].values,
                         '4':b['pt'].values,# + diff,
                         '5':b['pc'].values,# + diff,
                         '6':b['pn'].values,# + diff,
                         '7':b['mu'].values,
                         '8':b['crl'].values,
                         '9':b['chisqr'].values}
                         )

    if save is True:
        folder = p['rsr_path'] + '/' + pst
        if not os.path.exists(folder):
            os.makedirs(folder)
        data.to_csv(folder + '/' + product+'.'+pik, sep='\t', float_format='%.7f',
                    na_rep='nan', header=False, index=False)
        print('CREATED: ' + folder + '/' + product+'.'+pik)



@timing
@loop
def topik1m(pst, pik, from_process='pik1', from_product='MagLoResInco1', to_product='MagHiResInco1'):
    """ Inteprolate any pik file into 1m sampling
        NOTE: at this point, from_process is mandatory for @loop to work
    """
    p = icp.get.params()

    source = p['pik_path'].replace( p['process'], '') + from_process + '/' + pst + '/' + from_product + '.' + pik
    bxds = p['cmp_path'] + '/' + pst + '/'+ to_product
    sweep = '/'.join([p['sweep_path'], pst, 'sweeps'])

    test = icp.read.isfile(source) * icp.read.isfile(bxds) * icp.read.isfile(sweep, verbose=False)
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
    os.system('rm -f ' + LU + ' ' + P)
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
      a['pst'] = np.full(len(r['xo']), pst)
      a['pik'] = np.full(len(r['xo']), pik)
      a['xo'] = xo
      a['crl'] = r['crl']
      a['longitude'] = icp.get.longitude(pst)[xo]
      a['latitude'] = icp.get.latitude(pst)[xo]
      a['roll'] = icp.get.roll(pst)[xo]
      a['surface_range'] = icp.get.surface_range(pst)[xo]
      a['ice_thickness'] = icp.get.ice_thickness(pst)[xo]
    else:
      return

    if os.path.isfile(p['rsr_path'] + '/' + pst + '/' + product + '.' + pik + '.surface_coefficients'):
        b = icp.read.surface_coefficients(pst, pik, **kwargs)
        a['Rsc'] = b['Rsc']
        a['Rsn'] = b['Rsn']
    else:
        a['Rsc'] = [np.nan for i in xo]
        a['Rsn'] = [np.nan for i in xo]

    if os.path.isfile(p['rsr_path'] + '/' + pst + '/' + product + '.' + pik + '.surface_properties'):
        b = icp.read.surface_properties(pst, pik, **kwargs)
        a['sh'] = b['sh']
        a['eps'] = b['eps']
#        a['flag'] = r['flag']*b['flag']
    else:
        a['sh'] = [np.nan for i in xo]
        a['eps'] = [np.nan for i in xo]
#        a['flag'] = [np.nan for i in xo]

    if os.path.isfile(p['rsr_path'] + '/' + pst + '/' + product + '.' + pik + '.bed_coefficients'):
        b = icp.read.bed_coefficients(pst, pik, **kwargs)
        a['Rbc'] = b['Rbc']
        a['Rbn'] = b['Rbn']
    else:
        a['Rbc'] = [np.nan for i in xo]
        a['Rbn'] = [np.nan for i in xo]

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
def surface_coefficients(pst, pik, wb=15e6, gain=0, product='MagHiResInco1', save=True, **kwargs):
    """Surface coefficients (Reflectance and Scattering)
    """
    p = icp.get.params()
    a = icp.read.rsr(pst, pik, product='MagHiResInco1', **kwargs)
    if a is not None:
        h = icp.get.surface_range(pst)[a['xo'].astype(int)]
        Rsc, Rsn = invert.srf_coeff(Psc=a['pc'], Psn=a['pn'], h0=h, wb=15e6)

        out = {'0_Rsc':Rsc + gain, '1_Rsn':Rsn + gain}

        if save is True:
            p = icp.get.params()
            target = p['rsr_path'] + '/' + pst + '/' + product + '.' + pik + '.' + inspect.stack()[0][3]
            icp.do.save(out, target)



@loop
@timing
def surface_properties(pst, pik, wf=60e6, product='MagHiResInco1', save=True, **kwargs):
    """Return surface permittivity and RMS height
    """
    a = icp.read.rsr(pst, pik, **kwargs)
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
def bed_coefficients(pst, bed_pik, srf_pik=None, att_rate=0., wf=60e6, wb=15e6, product='MagHiResInco1', save=True, **kwargs):
    """Return bed reflance and scattering coefficients
    att_rate is 1-way attenuation rate of the ice in dB/km
    """
    p = icp.get.params()

    # Get srf_pik
    foo, pik_list = icp.get.pik(pst)
    srf_pik_list = [i for i in pik_list if 'srf' in i]

    if srf_pik:
        if srf_pik in srf_pik_list:
            pass
        else:
            print('IGNORED: No ' + srf_pik + ' for '+pst)
            return
    else: # Choose first pik in the list
        if srf_pik_list:
            srf_pik = srf_pik_list[0]
        else:
            print('IGNORED: No srf pik for '+pst)
            return

    s = icp.read.rsr(pst, srf_pik, **kwargs)
    s_prop = icp.read.surface_properties(pst, srf_pik, **kwargs)
    b = icp.read.rsr(pst, bed_pik, **kwargs)

    h0 = icp.get.surface_range(pst)[b['xo'].astype(int)]
    h1 = icp.get.ice_thickness(pst)[b['xo'].astype(int)]
    n1 = np.sqrt(s_prop['eps'])
    sh = s_prop['sh']
    Q1 = 2 * h1/1e3 * att_rate

    Rbc, Rbn = invert.bed_coeff(Psc=s['pc'], Psn=s['pn'],
                                Pbc=b['pc'], Pbn=b['pn'],
                                n1=n1, sh=sh, h0=h0, h1=h1, Q1=Q1,
                                wb=wb, wf=wf)

    out = {'0_Rbc':Rbc, '1_Rbn':Rbn}

    if save is True:
        p = icp.get.params()
        target = p['rsr_path'] + '/' + pst + '/' + product + '.' + bed_pik + '.' + inspect.stack()[0][3]
        icp.do.save(out, target)


