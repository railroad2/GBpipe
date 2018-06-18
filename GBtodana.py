import numpy as np

# spectrum of TOD
# 1/f noise fit

def tod_fft(data, clk=200e6, Nds=100000):
    Fs = clk / Nds # sampling frequency
    dt = 1/Fs
    N = data.size
    freq = np.fft.fftfreq(N, d=dt)

    dataf = np.fft.fft(data)
    
    return freq, dataf

if __name__=='__main__':
    test_todfft() 
