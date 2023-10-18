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
         linestyle = 'solid', linewidth = None, func = interp_linear, unique_marker='.'):
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
                p = plt.plot(x_space, equ(x_space), label = label, c = color, linestyle = linestyle, linewidth = linewidth, marker=unique_marker)
        else:
                p = plt.plot([], [], label = label, c = color, linestyle = linestyle, linewidth = linewidth, marker=unique_marker)

        plt.errorbar(x, y, xerr = xerr, yerr = yerr,
                        ms = marker_size, fmt = unique_marker,
                        c = p[-1].get_color())

        for i in exclude:
                plt.scatter(x[i], y[i], s = 60, marker = 'x', c = 'red')

        return pd.Series(x_clean), pd.Series(y_clean), equ


not_log_y1 = np.array([0.9925906573,0.981783889,0.9702959482,0.955498671,0.9386042172,0.9199664639,0.89991054,0.8787368754,0.8567108981,0.8342649048,0.7864084723,0.7982842488,0.6903632971,0.6468870705,0.606894041,0.4664580999,0.3739161264,0.3100183946,0.2638129759,0.2295238433,0.2027493198,0.1814937498,0.1641371804,0.1497575065,0.1056186798,0.05565599957,0.03710336106,0.0273750684,0.02137210158,0.01735252438,0.01641993018,0.01227514697,0.009952551212,0.008574197741,0.00788940987,0.00776859767,0.008326251145,0.009686663087,0.01372052189,0.02241504434,0.0151708743])
not_log_x1 = np.array([22.5,32.5,42.5,52.5,62.5,72.5,82.5,92.5,102.5,112.5,135,157.5,180,202.5,225,325,425,525,625,725,825,925,1025,1125,1600,3000,4400,5800,7200,8600,9000,11500,14000,16500,19000,21500,24000,26500,29000,31500,34000])

y1 = np.array([math.log(not_log_y1[i]) for i in range(len(not_log_y1))])
x1 = np.array([math.log(not_log_x1[i]) for i in range(len(not_log_x1))])

mu0 = 4 * math.pi * 10**(-7)
D = 45 * 10**(-3)
h = 1.5 * 10**(-3)
a = D / 2 - h

def ln_h1_h0(nu, sigma):
        delta = 1 / (sigma * mu0 * math.pi * nu)**0.5
        alpha = complex(1, 1) / delta
        # print(abs(np.cosh(alpha * h) + 0.5*a*alpha * np.sinh(h * alpha)))
        return np.log(1 / abs(np.cosh(alpha * h) + 0.5*a*alpha * np.sinh(h * alpha)))

x_th = np.linspace(0.5, 35000, 1000)

sigma_low = 4.72*10**7
sigma_mid = 5.89*10**7
sigma_high = 4.51*10**7
sigma_induct = 4.52*10**7

plt.figure(figsize = (15,8))

plt.yticks(np.arange(-6.5, 0.6, step=0.5))
# plt.xticks(np.arange(0, max(x1)*1.2, step=round(((max(x1) - min(x1))/len(x1)), 0)))

plt.ylim(min(y1)*1.2, 1)
plt.xlim(0, max(x1)*1.2)

plt.grid(color = 'black', linestyle = '--', linewidth = 0.5)

plt.title("$ln(\\frac{H_1}{H_0})(ln(ν))$")
plt.ylabel("$ln(\\frac{H_1}{H_0})$", fontsize=15)
plt.xlabel(f"ln(ν)", fontsize=15)

# print(ln_h1_h0(x_th, sigma_low))

plt.plot(np.log(x_th), ln_h1_h0(x_th, sigma_low))
plt.plot(np.log(x_th), ln_h1_h0(x_th, sigma_mid))
plt.plot(np.log(x_th), ln_h1_h0(x_th, sigma_high))
plt.plot(np.log(x_th), ln_h1_h0(x_th, sigma_induct))

plt.scatter(np.log(x_th)[::20], ln_h1_h0(x_th, sigma_low)[::20], label="$\sigma = {:.3g}$".format(sigma_low).replace("e+0", " \cdot 10^"), marker='^')
plt.scatter(np.log(x_th)[::20], ln_h1_h0(x_th, sigma_mid)[::20], label="$\sigma = {:.3g}$".format(sigma_mid).replace("e+0", " \cdot 10^"), marker='>')
plt.scatter(np.log(x_th)[::20], ln_h1_h0(x_th, sigma_high)[::20], label="$\sigma = {:.3g}$".format(sigma_high).replace("e+0", " \cdot 10^"), marker='*')
plt.scatter(np.log(x_th)[::20], ln_h1_h0(x_th, sigma_induct)[::20], label="$\sigma = {:.3g}$".format(sigma_induct).replace("e+0", " \cdot 10^"), marker='v')

x_clean, y_clean, equ = plot(x1, y1, x_min=0, func=None, label="Эксперимент")

plt.legend()

plt.savefig("gr6.png")
plt.show()
