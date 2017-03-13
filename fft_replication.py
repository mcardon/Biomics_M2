import random
from pylab import  conj, clf, plot, norm
from numpy import fft
import numpy as np


def get_template(M):
    W = [-0.5] * int(M/2) + [0.5] * int(M/2)
    return list(W * np.hanning(M))

def myfft(N=100000, M=1000, A=0.05, dispersion=0.01):
    """

    param N: length of the GC skew vector
    param M: length of the template
    param A: amplitude between positive and negative GC skew vector

    """
    N = int(N/2)
    x1 = [A+ (random.random()-0.5)*dispersion for i in range(N)]
    x2 = [-A+(random.random()-0.5)*dispersion for i in range(N)]
    x = x1 + x2
    x = x[int(N/4):] + x[0:int(N/4)]


    template =  get_template(M) + [0] * (N*2-M)
    template/=norm(template)

    c = abs(fft.ifft(
                    fft.fft(x) * conj(fft.fft(template))
                    )**2)/norm(x)/norm(template)

    # shift the SNR vector by the template length so that the peak is at the END of the template
    c = np.roll(c, M//2)

    return x, template, c*2./N


# Vary M to see if we need a full template or a short template
x, template, c = myfft(100000, 10000, 0.01); 
clf();
plot(x);
plot(template); plot(c*100000)

