## gvv_filter 7.1 for data
## Author : 
##--------Varun SM--------##

import soundfile as sf
import numpy as np
from scipy import signal


#Reading the soundfile
input_signal, fs = sf.read('Sound_Noise.wav')
sampl_freq = fs
order = 4
cutoff_freq = 4000 
x = np.pad(input_signal, (0,2**(int(np.log2(len(input_signal))) + 1) - len(input_signal)), 'constant', constant_values=0)
np.savetxt("../data/x.dat",np.array(x).reshape(len(x),1))
Wn = 2 * cutoff_freq / sampl_freq

#Passing butterworth filter
b, a = signal.butter(order, Wn, 'low')

n = int(len(a))
n = int(2 ** np.floor(np.log2(n)))
f = open("../data/a.dat", "w")
for i in range(n):
    f.write(str(a[i]) + "\n")
f.close()

n = int(len(b))
n = int(2 ** np.floor(np.log2(n)))
f = open("../data/b.dat", "w")
for i in range(n):
    f.write(str(b[i]) + "\n")
f.close()
z = np.array(x).reshape(len(x),1)

w = 2*np.pi*np.arange(len(z))/len(z)
z = np.exp(1j * w)
omega = np.linspace(0, 2*np.pi, len(input_signal))               # z = e^jw
num = np.polyval(b[::-1], z**(-1))
den = np.polyval(a[::-1], z**(-1))
H_k = num/den
f = open("../data/H_z.dat", "w")
for i in range(len(H_k)):
    f.write(str(H_k[i].real) + " " + str(H_k[i].imag) + "\n")
f.close()

