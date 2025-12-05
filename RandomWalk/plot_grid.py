import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

nodes_file = 'nodes.txt'
edges_file = 'edges.txt'
output_pdf = 'penrose_grid.pdf'


SMALL_SIZE = 14
MEDIUM_SIZE = 18
BIGGER_SIZE = 24
LEGEND_SIZE = 12
SIZE = 6

plt.rcParams['font.family'] = 'serif'
plt.rcParams["font.serif"] = ["Latin Modern Roman"]
plt.rcParams['mathtext.fontset'] = 'cm'

plt.rcParams['axes.titlepad'] = 10 
plt.rcParams['axes.labelpad'] = 10 
plt.rcParams['legend.fancybox'] = False
plt.rcParams['legend.edgecolor'] = "#000000"
plt.rcParams["legend.handlelength"] = 3
plt.rcParams["legend.framealpha"] = 1
plt.rcParams["legend.borderpad"] = 0.8

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=LEGEND_SIZE)   # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title



nodes = np.loadtxt(nodes_file)
edges = np.loadtxt(edges_file, dtype=int)
segments = nodes[edges - 1] 
fig, ax = plt.subplots(figsize=(12, 12))

lc = LineCollection(segments, colors='black', linewidths=0.5, alpha=0.7)

ax.add_collection(lc)


ax.set_title(r"Penrose Tiling Grid")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_aspect('equal') 
ax.autoscale_view()   

fig.savefig('penrose_grid.pdf', format='pdf', bbox_inches='tight')
plt.close(fig)
