#!env python3

import numpy as np
import wavio
import scipy.signal
import scipy.fftpack
import math
import sys
from scipy.fftpack import fft, ifft

# Generate the best linear approximation (BLA) of the system
# Note: The Fourier Transform IS the best least squares fit!
# Credit to Pascal Brunet:
# Nonlinear System Modeling and Identification of Loudspeakers
# https://repository.library.northeastern.edu/files/neu:336724/fulltext.pdf
def BLA(u, y, M, P, N, w):
    k = np.round(w/2/np.pi*N) + 1;
    K = len(k)
    Grz = np.zeros((K,M))
    Var_ng = np.zeros((K,M))
    realiz_u = np.reshape(u, ((*P+1* N), M))
    realiz_y = np.reshape(y, ((*P + 1 * N), M))

    for m in range(M):
        frames_u = np.reshape(realiz_u[N+1:, m], N, P)
        U = fft(frames_u)
        avgU = np.mean(U, axis=1)
        frames_y = np.reshape(realiz_y[N+1:,m], N, P)
        Y = fft(frames_y)
        avgY = np.mean(Y, axis=1)
        varY = np.var(Y, axis=1)
        Grz[:, m] = avgY[k] / avgU[k]
        Var_ng[:,m] = 1 / P*varY[k] / (np.abs(avgU[k]) ** 2)

    Gbla = np.mean(Grz,2)
    VarN = 1/M*np.mean(Var_ng, axis=1)
    VarBla = 1/M * np.var(Grz, axis=1)
    return Gbla, VarN, VarBla

