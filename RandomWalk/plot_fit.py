import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import ScalarFormatter

stats_file = 'Rstatistics.txt'
output_pdf = 'fit_plot.pdf'


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
plt.rcParams["legend.handlelength"] = 2.0
plt.rcParams["legend.framealpha"] = 1
plt.rcParams["legend.borderpad"] = 0.5

plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=MEDIUM_SIZE)
plt.rc('axes', labelsize=MEDIUM_SIZE)
plt.rc('xtick', labelsize=SMALL_SIZE)
plt.rc('ytick', labelsize=SMALL_SIZE)
plt.rc('legend', fontsize=LEGEND_SIZE)
plt.rc('figure', titlesize=BIGGER_SIZE)


data = np.loadtxt(stats_file)

n_steps = data[:, 0]
mean_s, std_s = data[:, 1], data[:, 2]
mean_nb, std_nb = data[:, 3], data[:, 4]
mean_oo, std_oo = data[:, 5], data[:, 6]


log_n = np.log10(n_steps)
log_R_s = np.log10(mean_s)
log_R_nb = np.log10(mean_nb)
log_R_oo = np.log10(mean_oo)

ln_10 = np.log(10)
err_s = std_s / mean_s / ln_10
err_nb = std_nb / mean_nb / ln_10
err_oo = std_oo / mean_oo / ln_10


w_s = 1.0 / (err_s**2)
w_nb = 1.0 / (err_nb**2)
w_oo = 1.0 / (err_oo**2)

def get_fit_params(x, y, w):
    fit, cov = np.polyfit(x, y, 1, w=w, cov=True)
    alpha = fit[0]
    err_alpha = np.sqrt(cov[0, 0])
    ln_c = fit[1]
    err_ln_c = np.sqrt(cov[1, 1])
    c = 10**ln_c
    err_c = c * err_ln_c * ln_10
    return alpha, err_alpha, c, err_c

alpha_s, err_alpha_s, c_s, err_c_s = get_fit_params(log_n, log_R_s, w_s)
alpha_nb, err_alpha_nb, c_nb, err_c_nb = get_fit_params(log_n, log_R_nb, w_nb)
alpha_oo, err_alpha_oo, c_oo, err_c_oo = get_fit_params(log_n, log_R_oo, w_oo)

print(f"  Simple RW     -> α = {alpha_s:.4f} ± {err_alpha_s:.4f}, c = {c_s:.3f} ± {err_c_s:.3f}")
print(f"  No-Back RW    -> α = {alpha_nb:.4f} ± {err_alpha_nb:.4f}, c = {c_nb:.3f} ± {err_c_nb:.3f}")
print(f"  Self-Avoiding -> α = {alpha_oo:.4f} ± {err_alpha_oo:.4f}, c = {c_oo:.3f} ± {err_c_oo:.3f}")


fit_x = np.geomspace(n_steps.min(), n_steps.max(), 300)
fit_y_s = c_s * (fit_x ** alpha_s)
fit_y_nb = c_nb * (fit_x ** alpha_nb)
fit_y_oo = c_oo * (fit_x ** alpha_oo)


fig, ax = plt.subplots(figsize=(10, 8))


ax.errorbar(n_steps, mean_s, yerr=std_s, fmt='o', capsize=5, color='red', label=f'Simple RW ($\\alpha = {alpha_s:.3f} \pm {err_alpha_s:.3f}$, $c = {c_s:.3f} \pm {err_c_s:.3f}$)')

ax.errorbar(n_steps, mean_nb, yerr=std_nb, fmt='s', capsize=5, color='blue', label=f'No-Back RW ($\\alpha = {alpha_nb:.3f} \pm {err_alpha_nb:.3f}$, $c = {c_nb:.3f} \pm {err_c_nb:.3f}$)')

ax.errorbar(n_steps, mean_oo, yerr=std_oo, fmt='^', capsize=5, color='green', label=f'Self-Avoiding ($\\alpha = {alpha_oo:.3f} \pm {err_alpha_oo:.3f}$, $c = {c_oo:.3f} \pm {err_c_oo:.3f}$)')


ax.plot(fit_x, fit_y_s, color='red', linestyle='-.', alpha=0.7)
ax.plot(fit_x, fit_y_nb, color='blue', linestyle='-.', alpha=0.7)
ax.plot(fit_x, fit_y_oo, color='green', linestyle='-.', alpha=0.7)

ax.set_xscale('log')
ax.set_yscale('log')


ax.set_title('Log-Log plot of R vs. n (Penrose Tiling)')
ax.set_xlabel('$n$ (Number of steps)')
ax.set_ylabel('$R$ (Mean distance)')


ax.legend(loc='best', frameon=False)

fig.savefig(output_pdf, format='pdf', bbox_inches='tight')
plt.close(fig)
