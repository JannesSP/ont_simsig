# Only executable on draco

import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy.optimize import curve_fit
from scipy.special import factorial
from scipy.stats import poisson

def fit_function(k, lamb):
    '''poisson function, parameter lamb is the fit parameter'''
    return poisson.pmf(k, lamb)

f = '/home/yi98suv/projects/modbuster/data/epinano/nanopolish/nanopolish_segmentation_bases_2.hdf5'
id_file = '/home/yi98suv/projects/modbuster/data/epinano/nanopolish/ids_nomod_rep1.ids'

ids = [l.strip() for l in open(id_file)]

h5 = h5py.File(f, 'r')

ms = []
md = []

for id in ids:
    if id in h5:
        ms.append(np.mean(np.diff(h5[id][:, 0])))
        md.append(np.median(np.diff(h5[id][:, 0])))

# ============================ MEANS ============================
# the bins should be of integer width, because poisson is an integer distribution
entries, bin_edges, patches = plt.hist(ms, bins=60, density=True, label='means')
axes = plt.gca()
y_min, y_max = axes.get_ylim()
plt.title('epinano_nomod_mean_segment_lengths_distribution.png')
plt.vlines(np.median(ms), ymin = y_min, ymax=y_max)
plt.text(np.median(ms), y_max, 'median: ' + str(np.median(ms)))

# calculate bin centres
bin_middles = 0.5 * (bin_edges[1:] + bin_edges[:-1])

# fit with curve_fit
parameters, cov_matrix = curve_fit(fit_function, bin_middles, entries)

# plot poisson-deviation with fitted parameter
x_plot = np.arange(0, 15)

plt.plot(
    x_plot,
    fit_function(x_plot, *parameters),
    marker='o', linestyle='',
    label='Poisson fit result',
)
plt.legend()
plt.savefig('epinano_nomod_mean_segment_lengths_distribution.png')
plt.close()


# ============================ MEDIAN ============================
# the bins should be of integer width, because poisson is an integer distribution
entries, bin_edges, patches = plt.hist(md, bins=60, density=True, label='medians')
axes = plt.gca()
y_min, y_max = axes.get_ylim()
plt.title('epinano_nomod_median_segment_lengths_distribution.png')
plt.vlines(np.median(md), ymin = y_min, ymax=y_max)
plt.text(np.median(md), y_max, 'median: ' + str(np.median(md)))

# calculate bin centres
bin_middles = 0.5 * (bin_edges[1:] + bin_edges[:-1])

# fit with curve_fit
parameters, cov_matrix = curve_fit(fit_function, bin_middles, entries)

# plot poisson-deviation with fitted parameter
x_plot = np.arange(0, 15)

plt.plot(
    x_plot,
    fit_function(x_plot, *parameters),
    marker='o', linestyle='',
    label='Poisson fit result',
)
plt.legend()
plt.savefig('epinano_nomod_median_segment_lengths_distribution.png')
plt.close()