import scipy.interpolate as inter
import numpy as np
import pylab as plt
from matplotlib.patches import FancyArrowPatch
import math

y1 = np.array([0.56,0.66,0.77,0.81,0.84,0.8954,0.78,0.81,0.72,0.87,0.79,0.53,0.82,0.803,0.622,0.7])
x1 = np.array([27.168,28.32,28.17,27.56,27.6,27.73,28.13,27.55,27.41,27.978,28.12,28.56,28.08,28.102,28.39,28.26])

y2 = np.array([0.269,0.25,0.275,0.222,0.29,0.253,0.298,0.23,0.3,0.229,0.272,0.2577,0.2278,0.205,0.2298,0.254])
x2 = np.array([19.378,19.07,19.68,18.83,19.41,19.96,19.54,18.844,19.57,18.82,19.08,18.98,18.79,18.69,18.78,18.91])


tetta_lvl = 1/np.sqrt(2)

ind = np.argsort(x1)
tck1,u1 = inter.splprep([x1[ind],y1[ind]], k = 3, s = 0.005)
out1 = inter.splev(u1,tck1)

ind = np.argsort(x2)
tck2,u2 = inter.splprep([x2[ind],y2[ind]], k = 3, s = 0.005)
out2 = inter.splev(u2,tck2)

plt.figure(figsize = (15,8))

# horizontal lines


plt.ylim(0, max(y1) + 0.5)
plt.xlim(min(x2) - 1, max(x1) + 1)

plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)

plt.title(f"Резонансные кривые в координатах U(f)")
plt.xlabel(f"f, кГц", fontsize=15)
plt.ylabel(f"U, В", fontsize=15)

plt.plot(out1[0], out1[1], label="При C₂", c = '#0C7EFA', linewidth=3, marker = 'o', mfc='#4690b3', ms=8)
plt.plot(out2[0], out2[1], label="При C₅", c = '#FA1D00', linewidth=3, marker = '^', mfc='#FA1D00', ms=8)

plt.legend()

plt.savefig("gr1.png")
plt.show()
