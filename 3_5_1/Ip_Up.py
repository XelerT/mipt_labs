import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import math
import scipy as sp
import scipy.optimize as op
from scipy.interpolate import UnivariateSpline
import copy

# Adding library directory to python search path
import sys
import os

import miptlib as mipt

# Create directory for generated .tex and .pdf
if not os.path.exists('gen'):
    os.mkdir('gen')


class consts:
    pi = 3.141592
    c  = 3 * 10**8 # m/s
    mu_0 = 4 * pi * 10**-7 # N * A^-2
    e = 1.6 * 10**-19 # C

    k = 1.38 * 10**-23 # J/K

class cgs:
    k = 1.38 * 10**-16 # erg/K
    e = 4.8 * 10**-10 # Franklin

class setup:
    dummy = 0

    d = 0.2 # mm
    l = 5.2 # mm
    S = d * l * consts.pi # mm^2

    U_lighting = 22.99 # V

    A1_div = 15 / 150 # mA / div
    err_A1 = 0.002 * 15 # mA
    m_i = 22 * 1.66 * 10**-27 # kg

    def U_real(x):
        return x

def sort_pairs(x, y, index):
    tup = list(zip(x, y))
    tup.sort(key = lambda x: x[index])

    untup = list(zip(*tup))

    return untup[0], untup[1]

def interp(x, y):
    return UnivariateSpline(x, y, s=100)

def deriv(x, y):
    return UnivariateSpline(x, y, s=100).derivative()

I_data = [2.42,3.01,3.35,3.7,3.04,4.3,4.68,4.97,4.74,4.38,4.04,3.64,3.05,2.55,2.17,1.92,1.82,1.64,1.39,1.1,0.88,0.72,0.51]
U_data = [22.56,20.87,20.14,19.99,19.66,19.53,19.53,19.35,19.5,19.47,19.68,19.96,20.67,22.17,23.24,24.2,24.68,25.48,27.01,30.17,33.46,34.39,35.19]

x, y = sort_pairs(setup.U_real(U_data), I_data, 0)
mipt.plot(x, y, func=interp)

max_der = deriv(x, y)(max(x))
print(max_der)
print("R_diff = ({:.3f} +- {:.3f}) Ом".format(1/max_der * 1000, 14))

plt.ylabel(r"$I_p$, мА")
plt.xlabel(r"$U_p$, В")
plt.grid(True)
#5 * 10^-3 = 14
# plt.legend()
plt.savefig("gen/Ip_Up.png")
    
