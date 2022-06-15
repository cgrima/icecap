# Calibration values

Calibration values for various seasons are reported here until it is hardcoded somewhere else in the hierarchy

The value is the *gain* keyword with which the epsilon value in the report is obtained.
Note that the Small Perturbation Method with the Small Angle Approximation must apply for epsilon to be correct (i.e., you should havem at least sh < 0.25 @ 60 MHz)

## ICP4

Calibration over blue ice (eps=3.15) at MIS

```python
In [1]: import icecap as icp                                                                                        
In [2]: icp.do.rsr('MIS/JKB2e/Y37a', 'srf_cyg', [46500,47999], air_loss=True, gain=-274.25).report()                               
[1] [ 11.52 s.] [ 33 eval.] [0.936]
Fit succeeded.
pt = -11.1 dB, pc = -11.1 dB, pn = -32.5 dB, pc-pn = 21.4 dB, 
SPM @ 60 MHz gives, eps = 3.147, sh = 3.50e-02 m
```

## ICP5

Calibration over blue ice (eps=3.15) at MIS

```python
In [1]: import icecap as icp
In [2]: icp.do.rsr('WLB/JKB2h/Y37b', 'srf_elg', [69815,70814], air_loss=True, gain=-274.20).report()                                                                    
[1] [  8.95 s.] [ 41 eval.] [0.947]
Fit succeeded.
pt = -11.1 dB, pc = -11.1 dB, pn = -34.3 dB, pc-pn = 23.2 dB, 
SPM @ 60 MHz gives, eps = 3.150, sh = 3.00e-02 m
```

## ICP6

Calibration over blue ice (eps=3.15) at MIS

```python
In [1]: import icecap as icp
Int [2]: icp.do.rsr('SMIS/MKB2l/X06a', 'srf_elg', [35000,36000], air_loss=True, gain=-270.23).report()
[1] [  7.73 s.] [ 43 eval.] [0.977]
Fit succeeded.
pt = -11.0 dB, pc = -11.2 dB, pn = -29.2 dB, pc-pn = 18.0 dB, 
SPM @ 60 MHz gives, eps = 3.149, sh = 5.00e-02 m
```

## ICP10_GCX

Calibration over blue ice patches (eps=3,15) at ELM1
```python
In [1]: import icecap as icp
In [2]: icp.do.rsr('ELM1/GCX0g/X44a', 'srf_lhb', [46750,47750], air_loss=True, gain=-259.17).report()
[1] [ 19.49 s.] [ 43 eval.] [0.920]
Fit succeeded.
pt = -11.3 dB, pc = -11.2 dB, pn = -26.1 dB, pc-pn = 14.8 dB, mu = 10.0 dB, 
SPM @ 60 MHz gives, eps = 3.148, sh = 7.49e-02 m
```

This gives the following permittivities over the calibration area:

![](https://github.com/cgrima/icecap/blob/master/figs/ICP10_GCX_calibration.png)

## GOG3

Pretty hard to find a good place for calibration with a stable aircraft. So, as for now, it is assumed to be the same as ICP5

## SRH1

Anja and I were pretty unsatisfied by the calibration we could do for SRH1. We ended up by using the surface snow density derived at CrA1 ice core (~ 137.8 kg/m3, i.e., permittivity = 1.246) since it is the only ice core with roughness smooth enough for our model to be reliable (< 0.25 m), but there is no way to measure the reliability of the derived calibration value, especially, I think the roll was pretty unstable. CrA1 and others are flown by DEV/JKB2t/Y75a. The figure below shows various RSR results for different calibration values.

![](https://github.com/cgrima/icecap/blob/master/figs/DEV_JKB2t_Y75a_calib.png)
