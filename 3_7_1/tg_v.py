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


def lin_plot(x, y, color, marker = '.', start=0, end=-1):
        coeffs = np.polyfit(x[start:end], y[start:end], 1)

        equ = lambda x: coeffs[0] * x + coeffs[1]
        
        plt.scatter(x, y, marker=marker, color = color, s=150)
        
        x_space = np.linspace(min(x)/5 - 100, max(x[start:end])*2, 100)    
        plt.plot(x_space, equ(x_space), color = color)

        return coeffs[0], coeffs[1] # k, b

y1 = np.array([0.9205648502,1,1.253960338,1.37638192,1.556030373,5.027339492,3.077683537,3.732050808,6.313751515,5.67128182,9.833802754,17.16933693])
x1 = np.array([135,157.5,180,202.5,225,325,425,525,625,725,825,925])

plt.figure(figsize = (15,8))

plt.grid(color = 'black', linestyle = '--', linewidth = 0.3)

plt.title("$tg(ψ)(ν)$")
plt.ylabel("$tg(ψ)$", fontsize=15)
plt.xlabel("$ν, Гц$", fontsize=15)

plt.yticks(np.arange(0, max(y1)*1.2, step=round(((max(y1) - min(y1))/len(y1)), 1)))
plt.xticks(np.arange(0, max(x1)*1.2, step=round(((max(x1) - min(x1))/len(x1)), 1)))

plt.ylim(0.7*min(y1), max(y1)*1.2)
plt.xlim(0, max(x1)*1.2)

k, b = lin_plot(x1, y1, 'black', end=5)
# plt.text(x=20, y=100, s="$δν = 1/T$", fontsize=17)

df, dat = least_sq(x1[0:5], y1[0:5])
latex = df.style.hide(level=0, axis=0).format(mnk_fmt).to_latex()
print(latex)

plt.savefig("gr2.png")
plt.show()

mu0 = 4 * math.pi * 10**(-7)
D = 45 * 10**(-3)
h = 1.5 * 10**(-3)
a = D / 2 - h
sigma = k / (math.pi * a * h * mu0)

print(f"k={k} b={b} sigma={sigma/10**(7)} *10^(7)")