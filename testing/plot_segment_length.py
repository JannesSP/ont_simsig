# Only executable on draco

import matplotlib.pyplot as plt
import numpy as np
import h5py
import os

path = os.path.join(os.path.dirname(__file__), '..', 'data', 'read_segment_length_distribution_epinano')

f = '/home/yi98suv/projects/modbuster/data/epinano/nanopolish/nanopolish_segmentation_bases_2.hdf5'
id_file = '/home/yi98suv/projects/modbuster/data/epinano/nanopolish/ids_nomod_rep1.ids'

ids = [l.strip() for l in open(id_file)]

h5 = h5py.File(f, 'r')

ms = []
md = []
mins = []

for id in ids:
    if id in h5:
        ms.append(np.mean(np.diff(h5[id][:, 0])))
        md.append(np.median(np.diff(h5[id][:, 0])))
        mins.append(np.diff(h5[id][:, 0]).min())

ms = np.array(ms)
md = np.array(md)
mins = np.array(mins)

# ============================ MEANS ============================
print('Plotting means')
# the bins should be of integer width, because poisson is an integer distribution
plt.title(f'EpiNano unmodified mean segment lengths distribution\nn={len(ms)}')
plt.hist(ms, bins=60, label='means', density=True)

# x = np.arange(min(ms), max(ms) + 1)

axes = plt.gca()
ymin, ymax = axes.get_ylim()
plt.vlines([np.mean(ms), np.median(ms)], ymin = ymin, ymax=ymax, color = 'red')
plt.text(np.mean(ms), ymax, 'mean: ' + str(np.mean(ms)))
plt.text(np.median(ms), ymax - .005, 'median: ' + str(np.median(ms)))

plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(path, 'epinano_nomod_mean_segment_lengths_distribution.png'))
plt.close()


# ============================ MEDIANS ============================
print('Plotting medians')
# the bins should be of integer width, because poisson is an integer distribution
plt.title(f'EpiNano unmodified median segment lengths distribution\nn={len(md)}')
plt.hist(md, bins=60, label='medians', density=True)

# x = np.arange(min(md), max(md) + 1)

axes = plt.gca()
ymin, ymax = axes.get_ylim()
plt.vlines([np.mean(md), np.median(md)], ymin = ymin, ymax=ymax, color = 'red')
plt.text(np.mean(md), ymax, 'mean: ' + str(np.mean(md)))
plt.text(np.median(md), ymax - .005, 'median: ' + str(np.median(md)))

plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(path, 'epinano_nomod_median_segment_lengths_distribution.png'))
plt.close()

# ============================ LOG MEANS ============================
print('Plotting logmeans')
logms = np.log(ms)
# the bins should be of integer width, because poisson is an integer distribution
plt.title(f'EpiNano unmodified log mean segment lengths distribution\nn={len(logms)}')
plt.hist(logms, bins=60, label='logmeans', density=True)

# x = np.arange(min(logms), max(logms) + 1)

axes = plt.gca()
ymin, ymax = axes.get_ylim()
plt.vlines([np.mean(logms), np.median(logms)], ymin = ymin, ymax=ymax, color = 'red')
plt.text(np.mean(logms), ymax, 'mean: ' + str(np.mean(logms)))
plt.text(np.median(logms), ymax - .005, 'median: ' + str(np.median(logms)))

plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(path, 'epinano_nomod_logmean_segment_lengths_distribution.png'))
plt.close()


# ============================ LOG MEDIANS ============================
print('Plotting logmedians')
logmd = np.log(md)
# the bins should be of integer width, because poisson is an integer distribution
plt.title(f'EpiNano unmodified log median segment lengths distribution\nn={len(logms)}')
plt.hist(logmd, bins=60, label='logmedians', density=True)

# x = np.arange(min(logmd), max(logmd) + 1)

axes = plt.gca()
ymin, ymax = axes.get_ylim()
plt.vlines([np.mean(logmd), np.median(logmd)], ymin = ymin, ymax=ymax, color = 'red')
plt.text(np.mean(logmd), ymax, 'mean: ' + str(np.mean(logmd)))
plt.text(np.median(logmd), ymax - .005, 'median: ' + str(np.median(logmd)))

plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(path, 'epinano_nomod_logmedian_segment_lengths_distribution.png'))
plt.close()

# ============================ MINIMUM ============================
print('Plotting minimum segment lengths')
# the bins should be of integer width, because poisson is an integer distribution
plt.title(f'EpiNano unmodified minimum segment lengths distribution\nn={len(mins)}min={min(mins)}')
plt.hist(mins, bins=60, label='minima', density=True)

# x = np.arange(min(logmd), max(logmd) + 1)

axes = plt.gca()
ymin, ymax = axes.get_ylim()
plt.vlines([np.mean(mins), np.median(mins)], ymin = ymin, ymax=ymax, color = 'red')
plt.text(np.mean(mins), ymax, 'mean: ' + str(np.mean(mins)))
plt.text(np.median(mins), ymax - .005, 'median: ' + str(np.median(mins)))

plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(path, 'epinano_nomod_min_segment_lengths_distribution.png'))
plt.close()

# ============================ DISTFIT ============================

# w = open('distfit_result.txt', 'w')

# dist = distfit(method='discrete')

# for data, string in zip([ms, md], ['means', 'medians']):
    
#     dist.fit_transform(data)
#     dist.plot()
#     plt.savefig(f'distfit_{string}.png')
#     plt.close()

#     w.write('string\n')
#     w.write(str(dist.model) + '\n')

# dist = distfit()

# for data, string in zip([logms, logmd], ['logmeans', 'logmedians']):
    
#     dist.fit_transform(data)
#     dist.plot()
#     plt.savefig(f'distfit_{string}.png')
#     plt.close()

#     w.write('string\n')
#     w.write(str(dist.model) + '\n')

# w.close()