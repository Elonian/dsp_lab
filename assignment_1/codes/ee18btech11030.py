## gvv_filter 7.1
## Author : 
##--------Varun SM--------##

import matplotlib.pyplot as plt
import numpy as np
import soundfile as sf
from scipy import signal


def FFT_vectorized(x):
    """A vectorized, non-recursive version of the Cooley-Tukey FFT"""
    x = np.asarray(x, dtype=float)
    N = x.shape[0]

    if np.log2(N) % 1 > 0:
        raise ValueError("size of x must be a power of 2")

    # N_min here is equivalent to the stopping condition above,
    # and should be a power of 2
    N_min = min(N, 32)

    # Perform an O[N^2] DFT on all length-N_min sub-problems at once
    n = np.arange(N_min)
    k = n[:, None]
    M = np.exp(-2j * np.pi * n * k / N_min)
    X = np.dot(M, x.reshape((N_min, -1)))

    # build-up each level of the recursive calculation all at once
    while X.shape[0] < N:
        X_even = X[:, :X.shape[1] / 2]
        X_odd = X[:, X.shape[1] / 2:]
        factor = np.exp(-1j * np.pi * np.arange(X.shape[0])
                        / X.shape[0])[:, None]
        X = np.vstack([X_even + factor * X_odd,
                       X_even - factor * X_odd])

    return X.ravel()

def new_filter(a, b, x):
    N = len(x)
    print(N)                                       # len of input               
    X_k = np.fft.fft(x)                              # fft of input
    omega = np.linspace(0, 2*np.pi, N)               # z = e^jw
    num = np.polyval(b[::-1], np.exp(1j*omega)**(-1)) 
    den = np.polyval(a[::-1], np.exp(1j*omega)**(-1))
    H_k = num/den                                    # H(e^jw) from coeff  
    Y_k = np.zeros(N) + 1j*np.zeros(N)
    for k in range(N):                               # Y(e^jw) = H(e^jw)*X(e^jw)
        Y_k[k] = X_k[k] * H_k[k]                         
    y = np.fft.ifft(Y_k)                             # ifft of Y(Z)
    return y,Y_k

input_signal,fs = sf.read('Sound_Noise.wav')         # Read the file
sampl_freq = fs                                      # fs of input
order = 4                                            # filter order
cutoff_freq = 4000.0                                 # cutoff freq
Wn = (2*cutoff_freq)/sampl_freq                      # digi freq
b, a = signal.butter(order, Wn, 'low')               # num, den coeff
output_signal = signal.filtfilt(b, a, input_signal)  # output from filter 
sf.write('Sound_Noise(Bfil).wav', output_signal, fs) # Create file and store
output_signal_freq = np.fft.fft(output_signal)

y,Y_k = new_filter(a, b, input_signal)
sf.write('Sound_Noise(nfil)_.wav', y.real, fs)


plt.title("Output from Butterworth filter(time)")    # Time domain plots
plt.grid()
plt.plot(output_signal,'y')
plt.savefig("../figs/ee18btech11030.eps")
plt.show()
plt.title("Output from own filter(time)")
plt.grid()
plt.plot(y,'k')
plt.savefig("../figs/ee18btech11030_1.eps")
plt.show()

plt.title("Output from Butterworth filter(freq)")    # Freq domain plots
plt.grid()
plt.plot(np.abs(np.fft.fftshift(output_signal_freq)), 'y')
plt.savefig("../figs/ee18btech11030_2.eps")
plt.show()
plt.title("Output from own filter(freq)")
plt.grid()
plt.plot(np.abs(np.fft.fftshift(Y_k)), 'k')
plt.savefig("../figs/ee18btech11030_3.eps")
plt.show()
