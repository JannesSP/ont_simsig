import os

import matplotlib.pyplot as plt
import pandas as pd
from numpy import diff

from Simulator import RNASimulator
from Writer import RNAWriter

# Loading model
kmer_model_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'template_median69pA.model')
df = pd.read_csv(kmer_model_file, sep='\t')
# kmers are stored from 3' to 5'
kmer_dict = {key : (mean, std) for key, mean, std in zip(df['kmer'], df['level_mean'], df['level_stdv'])}

print('Simulate RNA reads')
rna = RNASimulator(kmer_dict, length=20000, suffix='A'*50)
reference = rna.getReference()
num_of_signals = 12
signals = rna.drawRefSignals(num_of_signals)

writer = RNAWriter(reference, path=os.path.join('data','simulation'), dedicated_filename='test')
writer.writeReads(signals)

if not os.path.exists('plots'):
    os.makedirs('plots')

print('Plotting testdata')
fig, axs = plt.subplots(2, 3, figsize=(20,10))
for i, (signal, borders) in enumerate(signals[:3]):
    axs[0, i].plot(signal, alpha = 0.3)
    ylim = axs[0, i].get_ylim()
    xlim = axs[0, i].get_xlim()
    axs[0, i].text(xlim[0] + .05*(xlim[1]-xlim[0]), ylim[1] - .05*(ylim[1]-ylim[0]), '3\'', fontsize=13)
    axs[0, i].text(xlim[1] - .05*(xlim[1]-xlim[0]), ylim[1] - .05*(ylim[1]-ylim[0]), '5\'', fontsize=13)
    axs[1, i].hist(diff(borders), bins=100, alpha = 0.3, density = True)

axs[0, 0].set_ylabel('SIGNAL in pA')
axs[0, 3//2].set_xlabel('datapoints/time (~3 kHz)')

axs[1, 0].set_ylabel('FREQUENCY in density per read')
axs[1, 3//2].set_xlabel('segment length in datapoints')
    
fig.suptitle('SIMULATED RNA READS')
plt.tight_layout()
plt.savefig(os.path.join('plots','simulation.png'))