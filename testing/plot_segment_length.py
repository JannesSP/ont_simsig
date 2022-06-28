# Only executable on draco

import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy import stats
from distfit import distfit

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

ms = np.array(ms)
md = np.array(md)

# ============================ MEANS ============================

# the bins should be of integer width, because poisson is an integer distribution
plt.title('EpiNano unmodified mean segment lengths distribution')
plt.hist(ms, bins=60, label='means', density=True)

x = np.arange(min(ms), max(ms) + 1)

plt.plot(x, stats.poisson.pmf(x, np.mean(ms)), marker='o', label="poisson")

axes = plt.gca()
ymin, ymax = axes.get_ylim()
plt.vlines(np.median(ms), ymin = ymin, ymax=ymax, color = 'red')
plt.text(np.median(ms), ymax, 'median: ' + str(np.median(ms)))

plt.legend()
plt.savefig('epinano_nomod_mean_segment_lengths_distribution.png')
plt.close()


# ============================ MEDIANS ============================

# the bins should be of integer width, because poisson is an integer distribution
plt.title('EpiNano unmodified median segment lengths distribution')
plt.hist(md, bins=60, label='medians', density=True)

x = np.arange(min(md), max(md) + 1)

plt.plot(x, stats.poisson.pmf(x, np.mean(md)), marker='o', label="poisson")

axes = plt.gca()
ymin, ymax = axes.get_ylim()
plt.vlines(np.median(md), ymin = ymin, ymax=ymax, color = 'red')
plt.text(np.median(md), ymax, 'median: ' + str(np.median(md)))

plt.legend()
plt.savefig('epinano_nomod_median_segment_lengths_distribution.png')
plt.close()

# ============================ LOG MEANS ============================

logms = np.log(ms)
# the bins should be of integer width, because poisson is an integer distribution
plt.title('EpiNano unmodified log mean segment lengths distribution')
plt.hist(logms, bins=60, label='logmeans', density=True)

x = np.arange(min(logms), max(logms) + 1)

plt.plot(x, stats.poisson.pmf(x, np.mean(logms)), marker='o', label="poisson")

axes = plt.gca()
ymin, ymax = axes.get_ylim()
plt.vlines(np.median(logms), ymin = ymin, ymax=ymax, color = 'red')
plt.text(np.median(logms), ymax, 'median: ' + str(np.median(logms)))

plt.legend()
plt.savefig('epinano_nomod_logmean_segment_lengths_distribution.png')
plt.close()


# ============================ LOG MEDIANS ============================

logmd = np.log(md)
# the bins should be of integer width, because poisson is an integer distribution
plt.title('EpiNano unmodified log median segment lengths distribution')
plt.hist(logmd, bins=60, label='logmedians', density=True)

x = np.arange(min(logmd), max(logmd) + 1)

plt.plot(x, stats.poisson.pmf(x, np.mean(logmd)), marker='o', label="poisson")

axes = plt.gca()
ymin, ymax = axes.get_ylim()
plt.vlines(np.median(logmd), ymin = ymin, ymax=ymax, color = 'red')
plt.text(np.median(logmd), ymax, 'median: ' + str(np.median(logmd)))

plt.legend()
plt.savefig('epinano_nomod_logmedian_segment_lengths_distribution.png')
plt.close()

# ============================ DISTFIT ============================

w = open('distfit_result.txt', 'w')

dist = distfit(method='discrete')

for data, string in zip([ms, md], ['means', 'medians']):
    
    # normalize data counts to 1
    data = data / sum(data)
    dist.fit_transform(data)
    dist.plot()
    plt.savefig(f'distfit_{string}.png')
    plt.close()

    w.write('string\n')
    w.write(str(dist.model) + '\n')

dist = distfit()

for data, string in zip([logms, logmd], ['logmeans', 'logmedians']):
    
    # normalize data counts to 1
    data = data / sum(data)
    dist.fit_transform(data)
    dist.plot()
    plt.savefig(f'distfit_{string}.png')
    plt.close()

    w.write('string\n')
    w.write(str(dist.model) + '\n')

w.close()