# Only executable on draco

import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy import stats

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
plt.hist(ms, bins=60, density=True, label='means')
axes = plt.gca()
y_min, y_max = axes.get_ylim()
plt.title('epinano_nomod_mean_segment_lengths_distribution.png')
plt.vlines(np.median(ms), ymin = y_min, ymax=y_max)
plt.text(np.median(ms), y_max, 'median: ' + str(np.median(ms)))

# find minimum and maximum of xticks, so we know
# where we should compute theoretical distribution
xt = plt.xticks()[0]  
xmin, xmax = min(xt), max(xt)  
lnspc = np.linspace(xmin, xmax, len(ms))

plt.plot(lnspc, stats.poisson.pmf(lnspc, np.mean(ms)), label="poisson")
plt.plot(lnspc, stats.nbinom.pmf(lnspc, np.mean(ms), np.mean(ms)/2), label="nbinom")

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

# find minimum and maximum of xticks, so we know
# where we should compute theoretical distribution
xt = plt.xticks()[0]  
xmin, xmax = min(xt), max(xt)  
lnspc = np.linspace(xmin, xmax, len(md))

plt.plot(lnspc, stats.poisson.pmf(lnspc, np.mean(md)), label="poisson")
plt.plot(lnspc, stats.nbinom.pmf(lnspc, np.mean(md), np.mean(md)/2), label="nbinom")

plt.legend()
plt.savefig('epinano_nomod_median_segment_lengths_distribution.png')
plt.close()