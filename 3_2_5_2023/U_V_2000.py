import scipy.interpolate as inter
import numpy as np
import pylab as plt
from matplotlib.patches import FancyArrowPatch

y1 = np.array([0.412037037,0.537037037,0.7083333333,0.8796296296,0.9814814815,1,0.9722222222,0.9166666667,0.8564814815,0.8055555556])
x1 = np.array([0.625,0.703125,0.78125,0.859375,0.9375,1,1.09375,1.171875,1.25,1.328125])

tetta_lvl = 1/np.sqrt(2)

tck1,u1 = inter.splprep([x1,y1], k = 5, s = 1)
out1 = inter.splev(u1,tck1)


plt.figure(figsize = (15,8))

# horizontal lines

plt.axhline(y = tetta_lvl, xmin = 0.2, xmax = 0.9, color = 'black', linestyle = '--')     # ]
plt.axvline(x = 0.775, color = 'black', linestyle = '--')                                     #   to find 2ΔΩ
plt.axvline(x = 1.07, color = 'black', linestyle = '--')                                      # ] 

plt.arrow(x=0.775, y=tetta_lvl, dx=0.4, dy=0, head_width=0.04,head_length=0.01)
plt.arrow(x=1.05,  y=tetta_lvl, dx=-0.095, dy=0, head_width=0.04,head_length=0.01)
plt.text(x=1, y=tetta_lvl-0.1, s="2ΔΩ", fontsize=15)

print("2ΔΩ=", 1.07 - 0.775)
#

plt.ylim(0, max(y1) + 0.5)
plt.xlim(min(x1) - 0.05, max(x1) + 0.05)

plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)

plt.title(f"Резонансные кривые в координатах U/U₀(ν/ν₀)")
plt.xlabel(f"ν/ν₀", fontsize=15)
plt.ylabel(f"U/U₀", fontsize=15)

plt.plot(out1[0], out1[1], label="экспериментальная кривая", c = '#46a6b3', linewidth=3, marker = 'o', mfc='#4690b3', ms=8)

plt.savefig("U_V_2000.png")
plt.show()
