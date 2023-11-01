import scipy.interpolate as inter
import numpy as np
import pylab as plt
from matplotlib.patches import FancyArrowPatch
import math
import pandas as pd
import itertools

def interp_linear(x, y):
        coeffs = np.polyfit(x, y, 1)
        return lambda x: coeffs[0] * x + coeffs[1]

PLOT_MARKER = itertools.cycle(['.', 'v', '^', '<', '>', '*', 'o', '+', '1', '2', '3', '4'])
def plot(x, y, label = None, color = None, xerr = 0, yerr = 0,
         begin = 0, end = None, exclude = [],
         x_min = None, x_max = None, marker_size = 6,
         linestyle = 'solid', linewidth = None, func = interp_linear, unique_marker='.', mfc=None, ms=8):
        '''
        Creates plot with approximating line.

        :param x:                   x coordinates of points
        :param y:                   y coordinates of points
        :param label:               label of plot in legend
        :param color:               color of plot
        :param xerr:                x errors of points
        :param yerr:                y errors of points
        :param begin:               index of first point used for approximating line
        :param end:                 index of last point used for approximating line
        :param exclude:             indices of 'error' points
        :param x_min:               left end of approximating line
        :param x_max:               right end of approximating line
        :param func:                function for approximating, None -> no approximating line
        :param marker_size:         points size
        :param unique_marker:       True -> use internal unique marker, otherwise use unique_marker as marker itself

        :return x_clean, y_clean: pd.Series of values used for approximating line
        '''
        
        assert len(x) == len(y), "x and y must have same length"
        
        end = (len(x) - 1) if (end == None) else end
        
        x_clean = []
        y_clean = []
        for i in range(begin, end + 1):
                if i in exclude:
                        continue
                if np.isnan(x[i]) or np.isnan(y[i]):
                        continue
                x_clean.append(x[i])
                y_clean.append(y[i])
        
        x_min = min(x_clean) if (x_min == None) else x_min
        x_max = max(x_clean) if (x_max == None) else x_max

        x_space = np.linspace(x_min, x_max, 100)
        
        if unique_marker == True:
                unique_marker = next(PLOT_MARKER)
        
        equ = None
        p = None
        # At least two points and function for approximating.
        if (func != None and end - begin + 1 >= 2):
                equ = func(x_clean, y_clean)
                p = plt.plot(x_space, equ(x_space), label = label, c = color, linestyle = linestyle, linewidth = linewidth, marker=unique_marker, mfc=mfc, ms=ms)
        else:
                p = plt.plot([], [], label = label, c = color, linestyle = linestyle, linewidth = linewidth, marker=unique_marker, mfc=mfc, ms=ms)

        plt.errorbar(x, y, xerr = xerr, yerr = yerr,
                        ms = marker_size, fmt = unique_marker,
                        c = p[-1].get_color())

        for i in exclude:
                plt.scatter(x[i], y[i], s = 60, marker = 'x', c = 'red')

        return pd.Series(x_clean), pd.Series(y_clean), equ



y1 = np.array([0.5,0.5,0.7,0.9,0.65,0.5,1.5,1.5,1.6,1.5,1.285714286,1,0.67,0.5])
x1 = np.array([0.9350323044,0.9458004307,0.9669777459,0.9885139986,0.9623115578,0.9414931802,1.024766691,1.043072505,1.069633884,1.037329505,1.011844939,0.9996410625,0.9662598708,0.9411342426])

y2 = np.array([0.21,0.24,0.3,0.3333333333,0.3333333333,0.5,1,1.090909091,1.090909091,1.333333333,1.333333333,1.333333333])
x2 = np.array([0.91718107,0.9197530864,0.9393004115,0.9711934156,0.9763374486,0.9809670782,0.996399177,1.013888889,1.017489712,1.041666667,1.045781893,1.049382716])


tetta_lvl = 1/np.sqrt(2)

# ind = np.argsort(x1)
# tck1,u1 = inter.splprep([x1[ind],y1[ind]], k = 2, s = 1)
# out1 = inter.splev(u1,tck1)

# ind = np.argsort(x2)
# tck2,u2 = inter.splprep([x2[ind],y2[ind]], k = 3, s = 0.005)
# out2 = inter.splev(u2,tck2)

plt.figure(figsize = (15,8))

# horizontal lines

# plt.axhline(y = 3.237324758, xmin = -1, xmax = 100, color = 'black', linestyle = '--')     # ]
# plt.axvline(x = 0.981, color = 'black', linestyle = '--')                                     #   to find 2ΔΩ
# plt.axvline(x = 1.0215, color = 'black', linestyle = '--')                                      # ] 

# plt.arrow(x=1, y=tetta_lvl, dx=0.013, dy=0, head_width=0.02,head_length=0.005)
# plt.arrow(x=1,  y=tetta_lvl, dx=-0.017, dy=0, head_width=0.02,head_length=0.005)
# plt.text(x=0.995, y=tetta_lvl-0.1, s="2ΔΩ", fontsize=15)

# print("2ΔΩ=", 1.018-0.978)
#

plt.ylim(0, max(y1) + 0.5)
plt.xlim(min(x1) - 0.05, max(x1) + 0.05)

plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)

plt.title(f"Фазово­-частотные характеристики $\phi_U/\pi(f/f₀)$")
plt.xlabel(f"$f/f₀$", fontsize=15)
plt.ylabel(f"$\phi_U/\pi$", fontsize=15)


from scipy.interpolate import UnivariateSpline
def interp(x, y):
        return UnivariateSpline(x, y, s=0.3)

ind = np.argsort(x1)
x_clean, y_clean, equ = plot(x1[ind], y1[ind], func=interp, label="При C₂", color = '#0C7EFA', linewidth=3, unique_marker = 'o', mfc='#4690b3', ms=0)
ind = np.argsort(x2)
x_clean, y_clean, equ = plot(x2[ind], y2[ind], func=interp, label="При C₅", color = '#FA1D00', linewidth=3, unique_marker = '^', mfc='#FA1D00', ms=0)

plt.axhline(y = 3/4, xmin = 0, xmax = 10, color = 'black', linestyle = '--')
plt.axhline(y = 5/4, xmin = 0, xmax = 10, color = 'black', linestyle = '--')

plt.legend()

Q=1/(1.01-0.97)
print("Добротность C2=",Q)

Q=1/(1.03-0.99)
print("Добротность C5=",Q)


plt.savefig("gr3.png")
plt.show()
