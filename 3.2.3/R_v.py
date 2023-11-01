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



y1 = np.array([3.457394331,2.918733164,2.728211926,2.622722657,2.221999158])
x1 = np.array([32.09,27.77,23.182,21.157,15.969])

y2 = np.array([3.825980293,4.052595162,2.96934416,2.569873999,2.230167778])
x2 = np.array([32.079,27.813,23.143,21.18,16.08])


tetta_lvl = 1/np.sqrt(2)

plt.figure(figsize = (15,8))

plt.ylim(0, max(y2) + 0.5)
plt.xlim(15.969*0.6, max(x1) + 4)

plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)

plt.title(f"$R_L(υ_{{0n}})$")
plt.ylabel(f"$R_L$, Ом", fontsize=15)
plt.xlabel(f"$υ_{{0n}}$, кГц", fontsize=15)


from scipy.interpolate import UnivariateSpline
def interp(x, y):
        return UnivariateSpline(x, y, s=0.3)

ind = np.argsort(x1)
x_clean, y_clean, equ = plot(x1[ind], y1[ind], label="При $E_1$", color = '#0C7EFA', linewidth=3, unique_marker = 'o', mfc='#4690b3', ms=0)
ind = np.argsort(x2)
x_clean, y_clean, equ = plot(x2[ind], y2[ind], label="При $E_2$", color = '#FA1D00', linewidth=3, unique_marker = '^', mfc='#FA1D00', ms=0)

plt.legend()

plt.axhline(y = 2.707760099, xmin = -1, xmax = 100, color = '#0C7EFA', linestyle = '--')
plt.axhline(y = 3.0, xmin = -1, xmax = 100, color = '#FA1D00', linestyle = '-')

plt.savefig("gr4.png")
plt.show()
