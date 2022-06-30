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
signals = rna.drawRefSignals(20)

writer = RNAWriter(reference, path=os.path.join('data','simulation'), dedicated_filename='test')
writer.writeReads(signals)

if not os.path.exists('plots'):
    os.makedirs('plots')

figs, axs = plt.subplots(2, 1, figsize=(10,10))
for signal, borders in signals:
    axs[0].plot(signal, alpha = 0.5)
    axs[1].hist(diff(borders), bins=100, alpha = 0.5, density = True)

axs[0].set_xlabel('datapoints/time (~3 kHz)')
axs[0].set_ylabel('signal in pA')

axs[1].set_xlabel('segment length in datapoints')
axs[1].set_ylabel('Frequency in density per read')
plt.title('Simulated Reads')
plt.tight_layout()
plt.savefig(os.path.join('plots','simulation.png'))