## gvv_filter 7.1 for plotting
## Author : 
##--------Varun SM--------##
import numpy as np
import matplotlib.pyplot as plt
import soundfile as sf

x, fs = sf.read('Sound_Noise.wav')

output = np.loadtxt("../data/y.dat")
output = output[:, 0]

sf.write('Sound_Noise_c(fil).wav', output[0:len(x)], fs)
plt.figure(1)
plt.plot(output[0:len(x)], 'k')
plt.grid()
plt.title("Output signal using FFT in c")
plt.savefig("../figs/ee18btech11030_4.eps")
plt.show()
