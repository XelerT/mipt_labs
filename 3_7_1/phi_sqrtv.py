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

y1 = np.array([1.612133072,1.570796327,1.458596589,1.413716694,1.356596828,0.9817477042,1.099557429,1.047197551,0.9424777961,0.9599310886,0.8867398619,0,0.8435758051,0.7853981634,0.7853981634,0.8867398619,1.07763934,1.327052069,1.533396414,1.682996064,1.868705975,1.90739554,2.213394824,2.313740535,2.609548736,2.763437982,2.970853922,3.178992566,3.26560289,3.403392041,3.730641276])
x1 = np.array([11.61895004,12.5499004,13.41640786,14.23024947,15,18.02775638,20.61552813,22.91287847,25,26.92582404,28.72281323,30.41381265,32.01562119,0,33.54101966,40,54.77225575,66.33249581,76.15773106,84.85281374,92.73618495,94.86832981,107.2380529,118.3215957,128.4523258,137.8404875,146.628783,154.9193338,162.788206,170.2938637,177.4823935])

plt.figure(figsize = (15,8))

plt.grid(color = 'black', linestyle = '--', linewidth = 0.3)

plt.title("$(ψ - π/4)(\sqrt{ν})$")
plt.ylabel("$ψ - π/4$", fontsize=15)
plt.xlabel("$\sqrt{ν}$", fontsize=15)

plt.yticks(np.arange(-10, max(y1)*1.2, step=round(((max(y1) - min(y1))/len(y1)), 1)))
plt.xticks(np.arange(0, max(x1)*1.2, step=round(((max(x1) - min(x1))/len(x1)*2), 0)))

plt.ylim(0.7*min(y1), max(y1)*1.2)
plt.xlim(0, max(x1)*1.2)

k, b = lin_plot(x1, y1, 'black', start=15)
# plt.text(x=20, y=100, s="$δν = 1/T$", fontsize=17)

df, dat = least_sq(x1, y1)
latex = df.style.hide(level=0, axis=0).format(mnk_fmt).to_latex()
print(latex)

plt.savefig("gr3.png")
# plt.show()

mu0 = 4 * math.pi * 10**(-7)
D = 45 * 10**(-3)
h = 1.5 * 10**(-3)
a = D / 2 - h
sigma = k**2 / (math.pi * h**2 * mu0)

print(f"k={k} b={b} sigma={sigma/10**(7)} *10^(7)")