import scipy.interpolate as inter
import numpy as np
import pylab as plt

y1 = np.array([2.517546689,2.470045808,2.37755732,2.265967949,2.091546725,1.909585679,1.66077154,1.402406961,1.134491939,0.9505202733,0.7441804678,0.6273132211,0.5390972994,0.4146902303])
x1 = np.array([33.30088213,33.92920066,34.55751919,35.18583772,35.81415625,36.44247478,37.07079331,37.69911184,38.32743037,38.9557489,39.58406744,40.21238597,40.8407045,41.46902303])

y2 = np.array([0.6240459647,0.6715468456,0.7640353334,0.8756247044,1.050045929,1.232006975,1.480821113,1.739185693,2.007100715,2.19107238,2.397412186,2.514279433,2.602495354,2.726902423])
# x2 = np.array([33.30088213,33.92920066,34.55751919,35.18583772,35.81415625,36.44247478,37.07079331,37.69911184,38.32743037,38.9557489,39.58406744,40.21238597,40.8407045,41.46902303])

w0=37.69911184

tck1,u1 = inter.splprep([x1,y1], k = 3, s = 0.1)
out1 = inter.splev(u1,tck1)

tck2,u2 = inter.splprep([x1,y2], k = 3, s = 0.1)
out2 = inter.splev(u2,tck2)

plt.figure(figsize = (15,8))

# horizontal lines
plt.axhline(y = np.pi,     color = 'black', linestyle = '--')
plt.axhline(y = np.pi / 2, color = 'black', linestyle = '--')
plt.axhline(y = np.pi / 4, color = 'black', linestyle = '--')

pi_lines_values = np.array([np.pi, np.pi / 2, np.pi / 4])
pi_lines_names  = ['π', 'π/2', 'π/4']
# plt.yticks(pi_lines_values, pi_lines_names)
#

plt.axvline(x = w0, color = 'black', linestyle = '--')

plt.ylim(0, max(y1) + 1)
plt.xlim(min(x1) - 3, max(x1) + 3)

plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)

plt.title(f"Разность фаз от циклической частоты")
plt.xlabel(f"\u03C9, 1/с", fontsize=15)
plt.ylabel(f"\u0394\u03C6, рад", fontsize=15)

plt.plot(out1[0], out1[1], label="экспериментальная кривая", c = '#46a6b3', linewidth=3,                marker = 'o', mfc='#4690b3', ms=8)
plt.plot(out2[0], out2[1], label="отраженная кривая",        c = '#46b363', linestyle=':', linewidth=3, marker = 'x', mfc='#46b379', ms=8)

# plt.legend(loc ="upper right")
plt.legend(bbox_to_anchor =(0.8, 1.15), ncol=2, fontsize=13)

dW = 40.21238597 - 34.55751919
w0=37.69911184
Q=w0/dW
print("Добротность=",Q)

plt.savefig("gr5.png")
plt.show()