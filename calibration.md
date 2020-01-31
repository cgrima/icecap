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

## GOG3

Pretty hard to find a good place for calibration with a stable aircraft. So, is for now, t is assumed to be the same as ICP5
