import scipy.interpolate as inter
import numpy as np
import pylab as plt
from matplotlib.patches import FancyArrowPatch

y1 = np.array([0.227,0.118,0.072,0.042,0.038,0.035,0.031,0.027,0.027])
x1 = np.array([200,400,600,800,1000,1200,1400,1600,1800])

tetta_lvl = 1/np.sqrt(2)

tck1,u1 = inter.splprep([x1,y1], k = 5, s = 150)
out1 = inter.splev(u1,tck1)


plt.figure(figsize = (15,8))


plt.ylim(0, max(y1)*1.2)
plt.xlim(0, max(x1)*1.2)

plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)

plt.title(f"K(ν)")
plt.ylabel(f"K", fontsize=15)
plt.xlabel(f"ν, кГц", fontsize=15)

plt.plot(out1[0], out1[1], label="экспериментальная кривая", c = '#46a6b3', linewidth=3, marker = 'o', mfc='#4690b3', ms=8)


plt.savefig("gr4.png")
plt.show()
