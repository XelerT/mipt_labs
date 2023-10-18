import scipy.interpolate as inter
import numpy as np
import pylab as plt
from matplotlib.patches import FancyArrowPatch

y1 = np.array([10.658,7.65,5.6562,4.3054,3.8924,3.6456,3.3815,3.07,3.0374,3.0501,3.0745,3.1468,3.2349,4.0422,5.3903])
x1 = np.array([40,150,250,400,500,600,800,2000,4000,6000,7500,10000,12000,20000,25000])

plt.figure(figsize = (15,8))

plt.yticks(np.arange(-10, max(y1)*1.2, step=round(((max(y1) - min(y1))/len(y1)), 1)))
plt.xticks(np.arange(0, max(x1)*1.2, step=round(((max(x1) - min(x1))/len(x1)), 0)))

plt.ylim(0, max(y1)*1.2)
plt.xlim(0, max(x1)*1.2)

plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)

plt.title(f"L(ν)")
plt.ylabel(f"L, мГн", fontsize=15)
plt.xlabel(f"ν, Гц", fontsize=15)

plt.plot(x1, y1, label="экспериментальная кривая", c = '#46a6b3', linewidth=3, marker = 'o', mfc='#4690b3', ms=8)

plt.savefig("gr4.png")
# plt.show()
