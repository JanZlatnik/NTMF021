import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import brentq

data_file = 'results/percolation_data.txt'
output_pdf = 'results/percolation_curves.pdf'

SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 24
LEGEND_SIZE = 12

plt.rcParams['font.family'] = 'serif'
plt.rcParams["font.serif"] = ["Latin Modern Roman"]
plt.rcParams['mathtext.fontset'] = 'cm'

plt.rcParams['axes.titlepad'] = 10 
plt.rcParams['axes.labelpad'] = 10 
plt.rcParams['legend.fancybox'] = False
plt.rcParams['legend.edgecolor'] = "#000000"
plt.rcParams["legend.handlelength"] = 3
plt.rcParams["legend.framealpha"] = 0
plt.rcParams["legend.borderpad"] = 0.8

plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=MEDIUM_SIZE)
plt.rc('axes', labelsize=MEDIUM_SIZE)
plt.rc('xtick', labelsize=SMALL_SIZE)
plt.rc('ytick', labelsize=SMALL_SIZE)
plt.rc('legend', fontsize=LEGEND_SIZE)
plt.rc('figure', titlesize=BIGGER_SIZE)

data = np.loadtxt(data_file)
p = data[:, 0]
y7 = data[:, 1]
y9 = data[:, 2]
y11 = data[:, 3]
y13 = data[:, 4] 


mask = (p > 0.58) & (p < 0.59)
p_search = p[mask]
diff_search = np.abs(y11[mask] - y13[mask])
local_min_idx = np.argmin(diff_search)
idx = np.where(p == p_search[local_min_idx])[0][0]
xa = p[idx-1]
xb = p[idx+1]
y11a = y11[idx-1]
y11b = y11[idx+11]
y13a = y13[idx-1]
y13b = y13[idx+11]

m1 = (y11a - y11b) / (xa - xb)
m2 = (y13a - y13b) / (xa - xb)
pcrit = xa + (y13a - y11a) / (m1 - m2)

print(f"Critical percolation probability pc = {pcrit:.6f}")





fig, ax = plt.subplots(figsize=(10, 8))

ax.plot(p, y7, label=r'$N=7$', color='orange', linewidth=1.5)
ax.plot(p, y9, label=r'$N=9$', color='green', linewidth=1.5)
ax.plot(p, y11, label=r'$N=11$', color='blue', linewidth=1.5)
ax.plot(p, y13, label=r'$N=13$', color='red', linewidth=1.5)

ax.axvline(pcrit, color='black', linestyle='-.', linewidth=1, label=f'$p_c \\approx {pcrit:.4f}$', alpha=0.5)

ax.set_xlabel(r"Occupation probability")
ax.set_ylabel(r"Percolation probability")

ax.set_xlim(0, 1)
ax.set_ylim(-0.05, 1.05)

ax.legend()

fig.savefig(output_pdf, format='pdf', bbox_inches='tight')
plt.close(fig)

