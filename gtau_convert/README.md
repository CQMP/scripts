### gtau_convert
An executable and a python module that makes a Fourier transform from g(iw) to g(tau) with adding tails

#### Usage 
See example

```
import numpy as np
import pygtau_convert as pygt

gw_input = np.loadtxt("gw.dat", unpack=True)
gw = np.zeros([2,gw_input.shape[1]])*1j
gw[0] = gw_input[0] # Matsubara grid
gw[1] = gw_input[1] + 1j*gw_input[2] # g(iw)

gtau = pygt.convert(gw, 256, np.array([1,0,0])).transpose()
print gtau[0]
```
