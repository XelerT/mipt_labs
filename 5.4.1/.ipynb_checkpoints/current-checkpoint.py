import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statistics import mean

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
                    mean([513,  492,  516 ]),
                    mean([117,  112,  123 ]),
                    mean([513,  492,  516 ]),
                    mean([34,   26,   33  ]),
                    mean([2,    3,    3   ]),
                    mean([2,    3,    3   ]),
                    mean([5,    5,    3   ]),
                    mean([3,    4,    4   ]),
                    mean([2,    4,    4   ]),
                    mean([4,    3,    2   ]),
                    mean([3,   3,     3   ])])

I      =       np.array([0.05, 0.62, 1.45, 2.27, 3.14, 4.03, 4.95, 5.8, 7.01, 8.05, 9.10, 9.5, 9.96, 10,  10.2, 10.2,  10.15, 10.15,  10.05, 10, 9.95, 9.9])
pres   = 760 - np.array([735,   700, 650,  600,  550,  500,  450,  400, 350,  300,  250,  230, 210,  200, 180,  150,   120  , 100,    70,    50, 30,   0])

plot_1 = plt.figure(figsize=(12*0.9,10*0.9))
plt.grid(visible=True, linewidth=0.6)

plt.xlim(xmin=0, xmax=max(pres) * 1.1)
plt.ylim(ymin=0, ymax=max(I) * 2)
plt.xlabel('P, $\\text{мм.рт.ст}$', fontsize=20, rotation=0)

plt.ylabel('I, $\\text{пА}$', fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=16, size=10)
plt.grid(visible=True, linewidth=0.6)

plt.errorbar(
    pres, I,
    fmt='ro',
    linewidth=1,
    markersize=6,
    elinewidth=1,
    xerr=0.01,
    yerr=0.01,
    color='red',
    label = '',
    capsize=4
)


inc = np.polyfit(pres[0:12], I[0:12], deg=1)
print (mnk(pres[0:12], I[0:12]))
dec = np.polyfit(pres[14:len(pres)], I[14:len(pres)], deg=1)
print (mnk(pres[14:len(pres)], I[14:len(pres)]))


# increase line
x = np.linspace(0, 700, 10000)
y =  inc[0] * x + inc[1]
plt.plot(x, y,"k-.", label = 'y = ' + str(round(inc[0], 3)) + 'x' + str(round(inc[1], 3)))
#

# plato line
x = np.linspace(500, 800, 10000)
y =  dec[0] * x + dec[1]
plt.plot(x, y,"k--", label = 'y = ' + str(round(dec[0], 4)) + 'x' + '+' + str(round(dec[1], 3)))
#

plt.scatter(pres, I)
plt.legend(fontsize=18, markerscale = 2)


plt.title ('\nРис. 5. зависимость $I(P)$'  , size = 19)
plot_1.savefig('pics/I.png')
