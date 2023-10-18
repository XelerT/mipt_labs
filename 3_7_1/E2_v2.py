import scipy.interpolate as inter
import numpy as np
import pandas as pd
import pylab as plt
from matplotlib.patches import FancyArrowPatch
import math

mnk_fmt = {
        '$\overline{x}$' : "{:.2e}",
        '$\sigma_x^2$'   : "{:.2e}",
        '$\overline{y}$' : "{:.2f}",
        '$\sigma_y^2$'   : "{:.2e}",
        '$r_{xy}$'       : "{:.2e}",
        '$a$'            : "{:.2f}",
        '$\Delta a$'     : "{:.2f}",
        '$b$'            : "{:.2f}",
        '$\Delta b$'     : "{:.2f}",
}

def least_sq(x, y):
        sx = (x**2).mean() - (x.mean())**2
        sy = (y**2).mean() - (y.mean())**2
        rxy = (y*x).mean() - (y.mean() * x.mean())
        a = rxy / sx
        da = (1/(len(x) - 2) * (sy/sx - a**2))**(0.5)
        b = y.mean() - a * x.mean()
        db = da*(sx + (x.mean())**2)**(1/2)
        dat = pd.DataFrame({
                '$\overline{x}$' : [x.mean()],
                '$\sigma_x^2$'   : [sx],
                '$\overline{y}$' : [y.mean()],
                '$\sigma_y^2$'   : [sy],
                '$r_{xy}$'       : [rxy],
                '$a$'            : [a],
                '$\Delta a$'     : [da],
                '$b$'            : [b],
                '$\Delta b$'     : [db],
        })
        return dat, [x.mean(), sx, y.mean(), sy, rxy, a, da, b, db]


def lin_plot(x, y, color, marker = '.'):
        coeffs = np.polyfit(x, y, 1)

        equ = lambda x: coeffs[0] * x + coeffs[1]
        
        plt.scatter(x, y, marker=marker, color = color, s=150)
        
        x_space = np.linspace(min(x)/5 - 100, max(x)*2, 100)    
        plt.plot(x_space, equ(x_space), color = color)

        return coeffs[0], coeffs[1] # k, b

y1 = np.array([2868.679584,2932.179862,3002.022709,3095.724001,3208.170245,3339.476693,3489.98638,3660.199066,3850.825383,4060.826739])
x1 = np.array([506.25,1056.25,1806.25,2756.25,3906.25,5256.25,6806.25,8556.25,10506.25,12656.25])


plt.figure(figsize = (15,8))

plt.grid(color = 'black', linestyle = '--', linewidth = 0.3)

plt.title("$1/ξ^2(ν^2)$")
plt.ylabel("$1/ξ^2, 1/(В\cdotА\cdotГц)^2$", fontsize=15)
plt.xlabel("$ν^2, Гц^2$", fontsize=15)

plt.yticks(np.arange(0, max(y1)*1.2, step=((max(y1) - min(y1))/len(y1))))
plt.xticks(np.arange(0, max(x1)*1.2, step=((max(x1) - min(x1))/len(x1))))

plt.ylim(0.7*min(y1), max(y1)*1.2)
plt.xlim(0, max(x1)*1.2)

k, b = lin_plot(x1, y1, 'black')
# plt.text(x=20, y=100, s="$δν = 1/T$", fontsize=17)

df, dat = least_sq(x1, y1)
latex = df.style.hide(level=0, axis=0).format(mnk_fmt).to_latex()
print(latex)

plt.savefig("gr1.png")
# plt.show()

mu0 = 4 * math.pi * 10**(-7)
D = 45 * 10**(-3)
h = 1.5 * 10**(-3)
a = D / 2 - h
xi0 = math.sqrt(1 / b)
sigma = math.sqrt(k) / (a * h * mu0 * 2*math.pi / (2*xi0))

print(f"k={k} b={b} xi0={xi0} sigma={sigma/10**7} * 10^(7)")