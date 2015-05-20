import numpy as np
import pygtau_convert as pygt

wnmax = 512
beta = 10
gw = np.zeros([2,wnmax])*1j
gw[0] = (np.array(range(0,wnmax))*2 + 1)*np.pi/beta
gw[1] = 0.5 / (gw[0]*1j - 0.5) + 0.5 / (gw[0]*1j + 0.5)

gtau = pygt.convert(gw, 256, np.array([1,0,0])).transpose()
np.savetxt("gtau.dat", gtau.transpose())
