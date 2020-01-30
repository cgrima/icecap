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
