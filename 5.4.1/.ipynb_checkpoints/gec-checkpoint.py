import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statistics import mean
from operator import mul
from scipy.optimize import curve_fit

def mnk(x, y):
    k = ((x*y).mean() - x.mean() * y.mean())/((x**2).mean() - (x.mean())**2)
    sigma_k = np.sqrt(1/len(y)) * np.sqrt(((y**2).mean() - (y.mean())**2)/((x**2).mean() - (x.mean())**2) - k**2)

    b = y.mean() - k*x.mean()
    sigma_b = sigma_k * np.sqrt((x**2).mean() - (x.mean())**2)
    ret = [k, sigma_k, b, sigma_b]
    return ret

def get_delta(arr):
    res = []
    i = 0
    while i  < len(arr) - 1:
        res.append(arr[i] - arr[i + 1])
        i += 1
    return np.array(res)

v = np.array([
1226 / 85.3,
1269 / 93.0,
1352 / 85.6,
1415 / 85.4,
1337 / 85.4,
1401 / 85.5,
1313 / 85.6,
1242 / 85.2,
1227 / 85.1,
1184 / 85.2,
1228 / 85.2,
1170 / 85.4,
1072 / 86.3,
1033 / 85.4,
986 / 85.0 ,
919 / 85.2 ,
757 / 85.6 ,
521 / 87.0 ,
457 / 85.3 ,
302 / 86.1 ,
142 / 86.2 ,
 82 / 86.2 ,
 46 / 97.5 ,
 24 / 90.6
])

l  = 10 *np.array([
0    ,
0.25 ,
0.5  ,
1    ,
1.75 ,
2    ,
2.5  ,
3    ,
4    ,
5    ,
6    ,
6.25 ,
6.5  ,
7    ,
7.25 ,
7.5  ,
7.75 ,
8    ,
8.25 ,
8.5  ,
8.75 ,
9    ,
9.25 ,
9.5  ])

plot_1 = plt.figure(figsize=(12*0.9,10*0.9))
plt.grid(visible=True, linewidth=0.6)

plt.xlim(xmin=0, xmax=max(l) * 1.1)
plt.ylim(ymin=0, ymax=max(v) * 1.1)
plt.xlabel('l, $\\text{мм}$', fontsize=20, rotation=0)

plt.ylabel('$N$, $c^{-1}$', fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=16, size=10)
plt.grid(visible=True, linewidth=0.6)

plt.errorbar(
    l, v,
    fmt='D',
    linewidth=1,
    markersize=5,
    elinewidth=1,
    xerr=0.05,
    yerr=0.01,
    color='red',
    label = '',
    capsize=4
)



z = np.polyfit(l[13:23], v[13:23], deg=1)
print (mnk(l[13:23], v[13:23]))

polynom_1 = np.poly1d(z)
plt.plot(l, polynom_1(l),"k-.", linewidth=2, label = 'y = ' + str(round(z[0], 3)) + 'x' + '+' + str(round(z[1])))

z = np.polyfit(l, v, deg=7)
polynom_1 = np.poly1d(z)
plt.plot(l, polynom_1(l),"k--", linewidth=1, color = 'red', label = '$y = A_1x^7 + A_2x^6 + ... + c$')
x = np.linspace(56, 93, 10000)
plt.plot(x[0:9999], abs(get_delta(polynom_1(x))/get_delta(x) * 20),"k--", linewidth=1, color = 'blue', label = '$-20\\frac{d\\frac{N}{t}}{dl}$')

plt.legend(fontsize=18, markerscale = 2)

plt.title ('\nРис. 3. Зависимость $\\frac{N}{t}(l)$'  , size = 19)
plot_1.savefig('pics/gec.png')
