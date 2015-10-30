# !! THIS IS AN EXAMPLE. USE YOUR OWN VALUES !!

"""
Parameters for ICP4 data
Author: Cyril Grima <cyril.grima@gmail.com>
"""

"""
INSTRUMENT SPECIFICITIES
------------------------
"""

frq = 60e6 # central frequency [Hz]
bdw = 15e6 # bandwidth [Hz]
season = 'ICP4'

scale = 1/1000. # scale factor applied to the raw data
antgain = 9.2 # antenna gain [dB]
pt = 67. # transmitted power (67dBm = 5000W) [dB]
logain = -191.9 # low gain channel [dB]
higain = -237.9 # high gain channel [dB]


"""
PROCESS SPECIFICITIES
---------------------
"""
process = 'MagHiResInco1'
#process = 'MagHiResInco2to1'


"""
PATHS
-----
"""
code_path = '../../../../../code/xtra/'+season+'/RSR/pik1.1m.RADnh3'
rsr_path = '../../../../../targ/xtra/'+season+'/RSR/pik1.1m.RADnh3'
cmp_path = '../../../../../targ/xtra/'+season+'/CMP/pik1.1m.RADnh3'
pik_path = '../../../../../orig/xtra/'+season+'/PIK/pik1.1m.RADnh3'
#pik_path = '../../../../../targ/artl/icecap_sabrina_grl/rsr'
foc_path = '../../../../../targ/xtra/'+season+'/FOC/Best_Versions/S1_POS'
norm_path = '../../../../../targ/norm'


"""
ABSOLUTE CALIBRATION
--------------------
Calibration is done by adjusting the RSR roughness-corrected reflectance over a
reference zone of known permittivity. This leads to the abs_calib factor.

For ICP4, the reference is located over the ablation area of MIS.
we have used the frames [46500:47999] on MIS/JKB2e/Y37a and considered a
permittivity of ~3.15 (pure compact ice) for that particular location.
The center of this window has the following characteristics:

Frame  47000
Eastings     309009
Northings  -1276517
lat      -77.954811
lon      166.392057
rng     1029.077878
roll       0.324656

After calibration, report of the RSR over that window is:

[  5.06 s.] [ 35 eval.] [True] Tolerance seems to be too small.
a = 0.278, mu = 32.493, s = 0.016, pt = 0.078, crl = 0.977
pc = -11.1 dB, pn = -32.8 dB, pt = -11.1 dB, 
SPM @ 60 MHz gives, eps = 3.149, sh = 0.035 m
"""

abs_calib = 6.03 # Power [dB]

