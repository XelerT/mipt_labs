import scipy.interpolate as inter
import numpy as np
import pylab as plt
from matplotlib.patches import FancyArrowPatch

y1 = np.array([0.4117647059,0.4647058824,0.5294117647,0.6294117647,0.7352941176,0.8529411765,0.9470588235,1.0,0.9470588235,0.8823529412,0.8,0.7176470588,0.6411764706,0.5882352941])
x1 = np.array([0.8833333333,0.9,0.9166666667,0.9333333333,0.95,0.9666666667,0.9833333333,1.0,1.016666667,1.033333333,1.05,1.066666667,1.083333333,1.1])

tetta_lvl = 1/np.sqrt(2)

tck1,u1 = inter.splprep([x1,y1], k = 3, s = 0.005)
out1 = inter.splev(u1,tck1)


plt.figure(figsize = (15,8))

# horizontal lines

# plt.axhline(y = tetta_lvl, xmin = 0.346, xmax = 0.760, color = 'black', linestyle = '--')     # ]
# plt.axvline(x = 0.945, color = 'black', linestyle = '--')                                     #   to find 2ΔΩ
# plt.axvline(x = 1.07, color = 'black', linestyle = '--')                                      # ] 

plt.arrow(x=0.945, y=tetta_lvl, dx=0.115, dy=0, head_width=0.04,head_length=0.01)
plt.arrow(x=1.05,  y=tetta_lvl, dx=-0.095, dy=0, head_width=0.04,head_length=0.01)
plt.text(x=1, y=tetta_lvl-0.1, s="2ΔΩ", fontsize=15)

print("2ΔΩ=", 1.07 - 0.945)
#

plt.ylim(0, max(y1) + 0.5)
plt.xlim(min(x1) - 0.05, max(x1) + 0.05)

plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)

plt.title(f"Резонансные кривые в координатах U/U₀(ν/ν₀)")
plt.xlabel(f"ν/ν₀", fontsize=15)
plt.ylabel(f"U/U₀", fontsize=15)

plt.plot(out1[0], out1[1], label="экспериментальная кривая", c = '#46a6b3', linewidth=3, marker = 'o', mfc='#4690b3', ms=8)

w0=37.69911184
Q=w0/(1.05-0.945)
print("Добротность=",Q)

plt.savefig("gr3.png")
plt.show()
