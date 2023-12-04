import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import math
import scipy as sp
import scipy.optimize as op
from scipy.interpolate import UnivariateSpline
import copy
from scipy import interpolate
import sympy as sp

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
        return x * 10

def sort_pairs(x, y, index):
    tup = list(zip(x, y))
    tup.sort(key = lambda x: x[index])

    untup = list(zip(*tup))

    return untup[0], untup[1]


# def interp(x, y):
#     return interpolate.interp1d(x, y, 'cubic')
def interp(x, y):
    return UnivariateSpline(x, y, s=100)

arr_T_e = []
arr_n_e = []
arr_I_p = []

lists = {
    'I_p' : [],
    'kT_e' : [],
    'n_e' : [],
    'omega_p' : [],
    'r_de' : [],
    'r_d' : [],
    'N_d' : [],
    'alpha' : []
}

def plot_zond(index, color, left, right, begin, end, data_U_pos, data_I_pos, data_U_neg, data_I_neg):
    x = np.concatenate((data_U_pos, -data_U_neg))
    y = np.concatenate((data_I_pos, -data_I_neg))

    lists['I_p'].append(float(index))
    print(index)

    x, y = sort_pairs(x, y, 0)
    
    mipt.plot(x, y, func=interp, color=color, label=index + " мА", unique_marker=True, marker_size=4, linewidth=1.3)
    
    # Asymptotes
    data_asl = [x[0:left], y[0:left]]
    data_asr = [x[len(x) - right:len(x)], y[len(y) - right:len(y)]]

    asl = mipt.interp_linear(data_asl[0], data_asl[1])
    asr = mipt.interp_linear(data_asr[0], data_asr[1])
    
    mipt.plot(data_asl[0], data_asl[1],
              color=color, linestyle='dashed', x_min=min(x), x_max=0, linewidth = 0.8, unique_marker=' ')

    mipt.plot(data_asr[0], data_asr[1],
              color=color, linestyle='dashed', x_min=0, x_max=max(x), linewidth = 0.8, unique_marker=' ')

    I_n = np.average([-asl(0), asr(0)])
    print("\tI_n = {:.3} mkA".format(I_n))

    # Tangent
    begin += len(x)//2
    end   += len(x)//2

    data_tan = [x[begin:end], y[begin:end]]
    tan = mipt.interp_linear(data_tan[0], data_tan[1])

    dI_dU = (tan(1) - tan(0)) / 1
    print("\tdI/dU(0) = {:.3} mkA/V".format(dI_dU))

    # Intersection
    line_tan = sp.Line((0, tan(0)), (1, tan(1)))
    line_Il = sp.Line((0, asl(0)), (1, asl(0)))
    line_Ir = sp.Line((0, asr(0)), (1, asr(0)))

    inter_l = line_tan.intersection(line_Il)[0]
    inter_r = line_tan.intersection(line_Ir)[0]

    dU_l = float(inter_l.x)
    dU_r = float(inter_r.x)
    print("\tleft  dU = {:.3} V".format(dU_l))
    print("\tright dU = {:.3} V".format(dU_r))

    # Tangent
    mipt.plot(data_tan[0], data_tan[1], color=color, linestyle='dotted', unique_marker=' ',
              x_min=1.5 * float(inter_l.x), x_max=1.5 * float(inter_r.x))

    # Horizontal on I_n
    plt.axhline(y=asl(0), color=color, linewidth=0.8, linestyle='dotted')
    plt.axhline(y=asr(0), color=color, linewidth=0.8, linestyle='dotted')

    # Intersections
    plt.scatter(x=inter_r.x, y=inter_r.y, color='r', marker='x', s=20)
    plt.scatter(x=inter_l.x, y=inter_l.y, color='r', marker='x', s=20)

    kT_e1 = 0.5 * I_n / dI_dU
    print("\tkT_e (I_n) = {:.3} eV".format(kT_e1))
    kT_e2 = 0.5 * np.average([-dU_l, dU_r])
    print("\tkT_e (dU) = {:.3} eV".format(kT_e2))

    lists['kT_e'].append(kT_e1)

    T_e = kT_e1 * 11800
    print("\tT_e = {:.3} K".format(T_e))

    n_e = (I_n * 10**-6) / (0.4 * consts.e * setup.S * 10**-6) * np.sqrt(setup.m_i / (2 * consts.k * T_e))
    print("\tn_e = {:.3} 1/m^3".format(n_e))
    lists['n_e'].append(n_e / 100**3)

    omega_p = 5.6 * 10**4 * np.sqrt(n_e / 100**3)
    print("\tomega_p = {:.3} rad/s".format(omega_p))
    lists['omega_p'].append(omega_p)

    r_de = np.sqrt(cgs.k * T_e / (8 * consts.pi * n_e / 100**3)) / cgs.e #CGS
    print("\tr_de = {:.3} cm".format(r_de))
    lists['r_de'].append(r_de)

    T_i = 300 # K
    r_d = np.sqrt(cgs.k * T_i / (8 * consts.pi * n_e / 100**3)) / cgs.e #CGS
    print("\tr_d = {:.3} cm".format(r_d))
    lists['r_d'].append(r_d)

    N_d = 4/3 * consts.pi * r_d**3 * (n_e / 100**3) # n_i = n_e
    print("\tN_d = {:.3}".format(N_d))
    lists['N_d'].append(N_d)

    n = (2/760 * 10**5) / (consts.k * T_i)
    print("\tn = {:.3} 1/m^3".format(n))

    alpha = n_e / n
    print("\talpha = {:.3}".format(alpha))
    lists['alpha'].append(alpha)

    arr_T_e.append(T_e)
    arr_n_e.append(n_e)

    print("")

    
data_U_pos5 = np.array([24.96,22.1,19.02,16.09,13.02,10.06,8,6.04,4.03,2.04])
data_I_pos5 = np.array([101.26,106.74,104.29,99.76,89.94,76.56,64.19,50.17,33.33,15.23])
data_U_neg5 = np.array([24.96,22.04,19.06,16.03,13.02,10.06,8.02,6.12,4.07,2.05])
data_I_neg5 = np.array([88.56,95.14,92.61,88.46,79.01,65.5,53.2,39.4,21.81,2.93])

data_U_pos3 = np.array([24.94,22.05,19.04,16.06,13.01,10.02,8.03,6.08,4.05])
data_I_pos3 = np.array([74.4,73.69,71.5,68.21,62.8,54.43,46.94,37.9,26.3])
data_U_neg3 = np.array([24.98,22.1,19,16.08,13.07,10.09,8,6.03,4.05,2.03,0.55])
data_I_neg3 = np.array([65.3,64.84,62.7,59.81,54.72,46.48,38.59,29.3,18.19,5.62,-3.68])

data_U_pos15 = np.array([24.95,22.07,19.02,16.05,13.05,10,8.05,6.05,4.02,2.03])
data_I_pos15 = np.array([39.35,38.02,36.75,35.27,33.12,29.45,26.02,21.45,15.57,9.09])
data_U_neg15 = np.array([25.03,22.03,19.02,16.03,13.05,10.04,8.074,6.017,4.07,2.06])
data_I_neg15 = np.array([33.45,32.24,31.17,29.88,27.94,24.55,21.21,16.6,11.16,4.6])

plot_zond("1.5", 'Green', 4, 4, -4, 4, data_U_pos15, data_I_pos15, data_U_neg15, data_I_neg15)
plot_zond("3", 'Blue', 4, 4, -4, 4, data_U_pos3, data_I_pos3, data_U_neg3, data_I_neg3)
plot_zond("5", 'Orange', 5, 5, -4, 4, data_U_pos5, data_I_pos5, data_U_neg5, data_I_neg5)

arr_I_p.append(1.5)
arr_I_p.append(3)
arr_I_p.append(5)

plt.xlabel(r"$U$, В")
plt.ylabel(r"$I$, мкА")
plt.grid(True)

plt.legend()
plt.savefig('gen/Iz_Uz.png')

results = pd.DataFrame.from_dict(lists)
print(results)



lists = {
    'I_p' : [],
    'kT_e' : [],
    'n_e' : [],
    'omega_p' : [],
    'r_de' : [],
    'r_d' : [],
    'N_d' : [],
    'alpha' : []
}


fmt = {
    'I_p' :  [r'$I_p$, мА', '{:.1f}'],
    'kT_e' :  [r'$kT_e$, эВ', '{:.1f}'],
    'n_e' :  [r'$n_e \cdot 10^{9}$, $1 / \text{см}^3$', '{:.1f}', -9],
    'omega_p' :  [r'$\omega_p \cdot 10^{9}$, рад/с', '{:.1f}', -9],
    'r_de' :  [r'$r_{De}$, мкм', '{:.2f}', 4],
    'r_d' :  [r'$r_{D}$, мкм', '{:.2f}', 4],
    'N_d' :  [r'$N_{d}$', '{:.1f}'],
    'alpha' :  [r'$\alpha \cdot 10^{-6}$', '{:.3f}', 6],
}

# names_list = [fmt.get(item, item) for item in fmt]

# dat = pd.DataFrame({"-" : names_list, "ferrit" : ferrit_list, "fesi" : fesi_list, "feni" : feni_list})
# dat = dat.set_index('-').transpose()
# dat.insert(0, '--', ['Феррит (Fe-Ni-Zn)', 'Fe-Si', 'Fe-Ni'])

tab = mipt.table(results, fmt)
tab.rename().data
tab.to_latex('gen/tab_res.tex')
print(results)

fig = plt.figure()
ax = fig.add_subplot()

mipt.plot(arr_I_p, arr_T_e)

def yfmt(arg, dummy):
    return '${:.2e} $'.format(arg).replace('e+04', '')

# Convert to scientific notation
ax.yaxis.set_major_formatter(ticker.FuncFormatter(yfmt))
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)

plt.xlabel(r"$I_p$, $мА$")
plt.ylabel(r"$T_e$, $K ⋅ 10^4$")
plt.grid(True)

# plt.legend()
plt.savefig('gen/T_e_I_p.png', pad_inches=0.3, bbox_inches='tight')

fig = plt.figure()
ax = fig.add_subplot()

mipt.plot(arr_I_p, arr_n_e)

def yfmt(arg, dummy):
    return '${:.2e}$'.format(arg).replace('e+16', '')

# Convert to scientific notation
ax.yaxis.set_major_formatter(ticker.FuncFormatter(yfmt))
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)

plt.xlabel(r"$I_p$, $мА$")
plt.ylabel(r"$n_e$, $1/м^3 ⋅ 10^{16}$")
plt.grid(True)

# plt.legend()
plt.savefig('gen/n_e_I_p.png', pad_inches=0.3, bbox_inches='tight')