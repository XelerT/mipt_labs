import scipy.interpolate as inter
import numpy as np
import pandas as pd
import pylab as plt
from matplotlib.patches import FancyArrowPatch

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
    print(f"k={coeffs[0]}, b={coeffs[1]}")
    equ = lambda x: coeffs[0] * x + coeffs[1]
    
    plt.scatter(x, y, marker=marker, color = color, s=150)
    
    x_space = np.linspace(min(x)/5 - 100, max(x)*2, 100)    
    plt.plot(x_space, equ(x_space), color = color)

y1 = np.array([40,92,110,140,155,170,180,200])
x1 = np.array([65.0,90.4,110.1,141.5,154.8,167.1,178.5,199.4])


plt.figure(figsize = (15,8))

plt.grid(color = 'black', linestyle = '--', linewidth = 0.3)

plt.title("$Т_{϶ᴋᴄᴨ}(Т_{теор})$")
plt.xlabel("$T_{теор}, мкс$", fontsize=15)
plt.ylabel("$T_{эксп}, мкс$", fontsize=15)

plt.yticks(np.arange(-50, max(y1)*1.2, step=10))
plt.xticks(np.arange(0, max(x1)*1.2, step=10))

plt.ylim(-50, max(y1)*1.2)
plt.xlim(0, max(x1)*1.2)

lin_plot(x1, y1, 'black')
plt.text(x=20, y=100, s="$Т_{϶ᴋᴄᴨ} = 1.1T_{теор} - 20$", fontsize=17)

df, dat = least_sq(x1, y1)
latex = df.style.hide(level=0, axis=0).format(mnk_fmt).to_latex()
print(latex)


plt.savefig("gr1.png")
plt.show()
