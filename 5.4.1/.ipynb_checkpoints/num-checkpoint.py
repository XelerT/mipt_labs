import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statistics import mean
from operator import mul
from scipy.optimize import curve_fit


def get_delta(arr):
    res = []
    i = 0
    while i  < len(arr) - 1:
        print (arr[i] - arr[i + 1])
        res.append(arr[i] - arr[i + 1])
        i += 1
    return np.array(res)

def mnk(x, y):
    k = ((x*y).mean() - x.mean() * y.mean())/((x**2).mean() - (x.mean())**2)
    sigma_k = np.sqrt(1/len(y)) * np.sqrt(((y**2).mean() - (y.mean())**2)/((x**2).mean() - (x.mean())**2) - k**2)

    b = y.mean() - k*x.mean()
    sigma_b = sigma_k * np.sqrt((x**2).mean() - (x.mean())**2)
    ret = [k, sigma_k, b, sigma_b]
    return ret

p_a = 103000
t = 10
count   = np.array([mean([3710, 3729, 3843]),
                    mean([3275, 3318, 3242]),
                    mean([2554,	2428, 2439]),
                    mean([1427, 1418, 1339]),
                    mean([513, 492,   516]),
                    mean([117, 112,   123]),
                    mean([34, 26,     33]),
                    mean([2, 3,   3]),
                    mean([2, 3,   3]),
                    mean([5, 5,   3]),
                    mean([3, 4,   4]),
                    mean([2, 4,   4]),
                    mean([4, 3, 2]),
                    mean([4, 3, 2])])


pres   = (760 - np.array([735, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100]))

get_delta(count)

plot_1 = plt.figure(figsize=(12*0.9,10*0.9))
plt.grid(visible=True, linewidth=0.6)

plt.xlim(xmin=0, xmax=max(pres) * 1.1)
plt.ylim(ymin=0, ymax=max(count) * 1.1)
plt.xlabel('P, $\\text{мм.рт.ст}$', fontsize=20, rotation=0)

plt.ylabel('N, $num$', fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=16, size=10)
plt.grid(visible=True, linewidth=0.6)

plt.errorbar(
    pres, count,
    fmt='D',
    linewidth=1,
    markersize=7,
    elinewidth=1,
    xerr=0.01,
    yerr=0.01,
    color='red',
    label = '',
    capsize=4
)

z = np.polyfit(pres[2:5], count[2:5], deg=1)
polynom_1 = np.poly1d(z)
print(polynom_1(pres))
plt.plot(pres, polynom_1(pres),"k-.", linewidth=2, label = 'y = ' + str(round(z[0], 3)) + 'x' + '+' + str(round(z[1])))

z = np.polyfit(pres, count, deg=7)
polynom_1 = np.poly1d(z)
x = np.linspace(0, 320, 10000)
plt.plot(x, polynom_1(x),"k--", linewidth=1, color = 'red', label = '$y = A_1x^7 + A_2x^6 + ... + c$')
plt.plot(x[0:9999], abs(get_delta(polynom_1(x))/get_delta(x)) * 120,"k--", linewidth=1, color = 'blue', label = '$-120\\frac{dN}{dP}$')


plt.legend(fontsize=18, markerscale = 2)

plt.title ('\nРис. 4. зависимость $N(P)$'  , size = 19)
plot_1.savefig('pics/N.png')
