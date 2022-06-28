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
plt.title('EpiNano unmodified mean segment lengths distribution')
plt.hist(ms, bins=60, label='means', density=True)

x = np.arange(min(ms), max(ms) + 1)

plt.plot(x, stats.poisson.pmf(x, np.mean(ms)), 'go', label="poisson", color = 'orange')
plt.plot(x, stats.nbinom.pmf(x, np.mean(ms), np.mean(ms)/2), 'go', label="nbinom", color = 'red')

axes = plt.gca()
ymin, ymax = axes.get_ylim()
plt.vlines(np.median(ms), ymin = ymin, ymax=ymax, color = 'red')
plt.text(np.median(ms), ymax, 'median: ' + str(np.median(ms)))

plt.legend()
plt.savefig('epinano_nomod_mean_segment_lengths_distribution.png')
plt.close()


# ============================ MEDIAN ============================
# the bins should be of integer width, because poisson is an integer distribution
plt.title('EpiNano unmodified median segment lengths distribution')
plt.hist(md, bins=60, label='medians', density=True)

x = np.arange(min(md), max(md) + 1)

plt.plot(x, stats.poisson.pmf(x, np.mean(md)), 'go', label="poisson", color = 'orange')
plt.plot(x, stats.nbinom.pmf(x, np.mean(md), np.mean(md)/2), 'go', label="nbinom", color = 'red')

axes = plt.gca()
ymin, ymax = axes.get_ylim()
plt.vlines(np.median(md), ymin = ymin, ymax=ymax, color = 'red')
plt.text(np.median(md), ymax, 'median: ' + str(np.median(md)))

plt.legend()
plt.savefig('epinano_nomod_median_segment_lengths_distribution.png')
plt.close()