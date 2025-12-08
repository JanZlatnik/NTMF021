import numpy as np
import matplotlib.pyplot as plt

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

data = np.loadtxt("results/rng_test.txt")

mean_val = np.mean(data)
std_val = np.std(data)
expected_mean = 0.5
expected_std = np.sqrt(1/12)

print(f"Number of samples: {len(data)}")
print(f"Mean:              {mean_val:.6f} (Expected: {expected_mean})")
print(f"Std Dev:           {std_val:.6f}  (Expected: {expected_std:.6f})")

fig, ax = plt.subplots(figsize=(10, 6))

ax.hist(data, bins=100, density=True, color='#87CEEB', edgecolor='black', alpha=0.7, label='Generated Data')
ax.plot([0, 1], [1, 1], linewidth=2, color='red', linestyle='--', label='Theoretical Uniform')

ax.set_title(f"RNG Quality Test ({len(data)} samples)")
ax.set_xlabel("Value")
ax.set_ylabel("Probability Density")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.2)

ax.legend(loc='best', ncol=2)

fig.savefig("results/rng_histogram.pdf", format='pdf', bbox_inches='tight')
plt.close(fig)