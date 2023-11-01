import scipy.interpolate as inter
import numpy as np
import pylab as plt
from matplotlib.patches import FancyArrowPatch
import math

y1 = np.array([0.6174200662,0.7276736494,0.848952591,0.8930540243,0.9261300992,0.9872105843,0.8599779493,
0.8930540243,
0.7938257993,
0.9592061742,
0.8710033076,
0.5843439912,
0.9040793826,
0.8853362734,
0.6857772878,
0.7717750827])
x1 = np.array([0.9751615219,1.016511127,1.011127064,0.9892318737,0.9906676238,0.9953338119,1.009691314,0.9888729361,0.9838478105,1.004235463,1.009332376,1.025125628,1.007896626,1.008686289,1.01902369,1.014357502])

y2 = np.array([0.7219538379,0.6709608159,0.7380568975,0.5958132045,0.7783145464,0.6790123457,0.7997852925,0.6172839506,0.8051529791,0.6146001074,0.7300053677,0.691626409,0.6113794954,0.550187869,0.616747182,0.6816961889])
x2 = np.array([0.9968106996,0.9809670782,1.012345679,0.9686213992,0.9984567901,1.026748971,1.005144033,0.9693415638,1.006687243,0.9681069959,0.9814814815,0.9763374486,0.966563786,0.9614197531,0.9660493827,0.9727366255])


tetta_lvl = 1/np.sqrt(2)

ind = np.argsort(x1)
tck1,u1 = inter.splprep([x1[ind],y1[ind]], k = 3, s = 0.005)
out1 = inter.splev(u1,tck1)

ind = np.argsort(x2)
tck2,u2 = inter.splprep([x2[ind],y2[ind]], k = 3, s = 0.005)
out2 = inter.splev(u2,tck2)

plt.figure(figsize = (15,8))

# horizontal lines

# plt.axhline(y = tetta_lvl, xmin = 0.346, xmax = 0.760, color = 'black', linestyle = '--')     # ]
# plt.axvline(x = 0.981, color = 'black', linestyle = '--')                                     #   to find 2ΔΩ
# plt.axvline(x = 1.0215, color = 'black', linestyle = '--')                                      # ] 

plt.arrow(x=1, y=tetta_lvl, dx=0.013, dy=0, head_width=0.02,head_length=0.005)
plt.arrow(x=1,  y=tetta_lvl, dx=-0.017, dy=0, head_width=0.02,head_length=0.005)
plt.text(x=0.995, y=tetta_lvl-0.1, s="2ΔΩ", fontsize=15)

print("2ΔΩ=", 1.018-0.978)
#

plt.ylim(0, max(y1) + 0.5)
plt.xlim(min(x1) - 0.05, max(x1) + 0.05)

plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)

plt.title(f"Резонансные кривые в координатах U/U₀(f/f₀)")
plt.xlabel(f"f/f₀", fontsize=15)
plt.ylabel(f"U/U₀", fontsize=15)

plt.plot(out1[0], out1[1], label="При C₂", c = '#0C7EFA', linewidth=3, marker = 'o', mfc='#4690b3', ms=8)
plt.plot(out2[0], out2[1], label="При C₅", c = '#FA1D00', linewidth=3, marker = '^', mfc='#FA1D00', ms=8)

plt.legend()

f0=27.86
Q=1/(1.018-0.978)
print("Добротность C2=",Q)

Q=1/(1.0215-0.981)
print("Добротность C5=",Q)


plt.savefig("gr2.png")
plt.show()
