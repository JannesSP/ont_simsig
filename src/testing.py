import os

import matplotlib.pyplot as plt
import pandas as pd
from numpy import diff

from Simulator import RNASimulator
from Writer import RNAWriter

# =============== Loading model ===============
kmer_model_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'template_median69pA.model')
df = pd.read_csv(kmer_model_file, sep='\t')
# kmers are stored from 3' to 5'
kmer_dict = {key : (mean, std) for key, mean, std in zip(df['kmer'], df['level_mean'], df['level_stdv'])}

# =============== Simulate ===============

segment_length = 50
for stdev_scale in [0, 0.2, 0.4, 0.6, 0.8, 1.0]:

    print(f'LOOP Segment length: {segment_length} & stdev scale: {stdev_scale}')

    rna = RNASimulator(kmer_dict, length=2000, suffix='A'*50, seed = 0, stdev_scale=stdev_scale, set_segment_length=segment_length)
    reference = rna.getReference()
    num_of_signals = 4000

    print('Simulate RNA reads')
    signals = rna.drawRefSignals(num_of_signals)
    # signals = rna.drawReadSignals(num_of_signals, min_len=1000, max_len=2000)

    # =============== Writing fast5 ===============

    path = os.path.join(os.path.dirname(__file__), '..', 'data')
    writer = RNAWriter(reference, path=path, tag = 'simulation')
    writer.writeReads(signals)

    # =============== Plotting first 3 reads ===============
    print('Plotting testdata')

    fig, axs = plt.subplots(2, 3, figsize=(20,10))
    for i, (signal, borders) in enumerate(signals[:3]):
        axs[0, i].plot(signal)
        ylim = axs[0, i].get_ylim()
        xlim = axs[0, i].get_xlim()
        axs[0, i].text(xlim[0] + .05*(xlim[1]-xlim[0]), ylim[1] - .05*(ylim[1]-ylim[0]), '3\'', fontsize=13)
        axs[0, i].text(xlim[1] - .05*(xlim[1]-xlim[0]), ylim[1] - .05*(ylim[1]-ylim[0]), '5\'', fontsize=13)
        axs[0, i].text(xlim[0] + .4*(xlim[1]-xlim[0]), ylim[1] - .05*(ylim[1]-ylim[0]), f'Readlength: {len(borders) - 1}', fontsize=13)
        axs[1, i].hist(diff(borders), bins=100, density = True)

    axs[0, 0].set_ylabel('SIGNAL in pA')
    axs[0, 3//2].set_xlabel('DATAPOINTS (~3 kHz)')

    axs[1, 0].set_ylabel('FREQUENCY in density per read')
    axs[1, 3//2].set_xlabel('SEGMENT-LENGTH in datapoints')
        
    fig.suptitle(f'SIMULATED RNA READS\nReference length: {rna.getRefLength()}')
    plt.tight_layout()
    plt.savefig(os.path.join(writer.getFilename() + 'simulation.png'))
    plt.close()

    # =============== Plotting whole read lengths and segmentation distribution ===============
    lengths = [len(read[0]) for read in signals]
    plt.hist(lengths, bins = 100)
    plt.title('Signal lengths distribution')
    plt.xlabel('Signal length')
    plt.ylabel('Frequency')
    plt.tight_layout()
    plt.savefig(os.path.join(writer.getFilename() + 'signal_lengths.png'))
    plt.close()