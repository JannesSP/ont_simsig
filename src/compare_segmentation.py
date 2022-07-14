import os
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO


def parse() -> Namespace:

    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
        ) 
    
    parser.add_argument('sim_fast5', type = str, help='FAST5 file with simulated reads')
    parser.add_argument('nanopolish_segmentation_fast5', type = str, help='FAST5 file with nanopolish segmentation')
    parser.add_argument('fastq', type = str, help = 'Basecalled reads')
    parser.add_argument('outpath', type = str, help = 'Path to output files')

    return parser.parse_args()

def compare(sim5 : h5py.File, nano5 : h5py.File, fastqs : dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    
    # store key features for each read describing the segmentation and basecalling quality
    read_statistics = pd.DataFrame(columns=['segmentation type', 'read', 'data Available', 'read segments mean', 'read segments stdev', 'basecall segmentation ratio', 'mean error', 'stdev error', 'max error', 'min error'])
    read_statistics = read_statistics.astype(dtype = {'segmentation type':'str', 'read':'str', 'data Available':'bool', 'read segments mean':'float', 'read segments stdev':'float', 'basecall segmentation ratio':'float', 'mean error':'float', 'stdev error':'float', 'max error':'float', 'min error':'float'})

    # analyse each individual segment
    segment_statistics = pd.DataFrame(columns=['segmentation type', 'read', 'segment mean', 'segment stdev', 'segment length'])
    segment_statistics = segment_statistics.astype(dtype = {'segmentation type':'str', 'read':'str', 'segment mean':'float', 'segment stdev':'float', 'segment length':'int'})

    noSeg = 0

    for i, read in enumerate(sim5):

        if not 'read' in read:
            continue

        print(f'Analysing segmentation of read {read} - {i}/{len(sim5) - 1} ...    ', end='\r')

        # key in sim5, ONT FAST5 format
        readid = read.split('_')[-1]

        simBorders = sim5[read]['Raw/Borders'][:]
        nSimSegments = sim5[read]['Raw'].attrs['num_segments']
        r_offset = sim5[read]['channel_id'].attrs['offset']
        r_range = sim5[read]['channel_id'].attrs['range']
        r_digit =sim5[read]['channel_id'].attrs['digitisation']
        signal = (sim5[read]['Raw/Signal'][:] + r_offset) * r_range / r_digit
        assert (nSimSegments + 1) == len(simBorders)

        if readid in nano5:
            nanoBorders = nano5[readid][:, 0]
            e = getSegmentErrors(nanoBorders, simBorders)
        # Nanopolish has no segmentation for this read
        # maybe read could not be mapped with minimap2
        else:
            nanoBorders = [np.nan]
            e = [np.nan]
            noSeg += 1
    
        read : str = str(fastqs[readid].seq)
        # phred : list[int] = fastqs[readid].letter_annotations['phred_quality']

        sim_entry = pd.DataFrame({
            'segmentation type':['simulation'],
            'read':[readid],
            'data Available':[True],
            'read segments mean':[np.mean([np.mean(signal[simBorders[i]:simBorders[i+1]]) for i in range(nSimSegments)])],
            'read segments stdev':[np.mean([np.std(signal[simBorders[i]:simBorders[i+1]]) for i in range(nSimSegments)])],
            'basecall segmentation ratio':[len(read)/nSimSegments],
            'mean error':[np.nan],
            'stdev error':[np.nan],
            'max error':[np.nan],
            'min error':[np.nan]
        })

        nano_entry = pd.DataFrame({
            'segmentation type':['nanopolish'],
            'read':[readid],
            'data Available':[True if len(nanoBorders) > 1 else False],
            'read segments mean':[np.mean([np.mean(signal[nanoBorders[i]:nanoBorders[i+1]]) for i in range(len(nanoBorders) - 1)]) if len(nanoBorders) > 1 else np.nan],
            'read segments stdev':[np.mean([np.std(signal[nanoBorders[i]:nanoBorders[i+1]]) for i in range(len(nanoBorders) - 1)]) if len(nanoBorders) > 1 else np.nan],
            'basecall segmentation ratio':[len(read)/(len(nanoBorders) - 1) if len(nanoBorders) > 1 else np.Infinity],
            'mean error':[np.mean(e)],
            'stdev error':[np.std(e)],
            'max error':[np.max(e)],
            'min error':[np.min(e)]
        })

        read_statistics = pd.concat([read_statistics, sim_entry, nano_entry], ignore_index=True)

        if readid in nano5:

            for sim_i in np.random.choice(nSimSegments, size = 100, replace = False): #range(nSimSegments):

                base_sim_entry = pd.DataFrame({
                    'segmentation type':['simulation'],
                    'read':[readid],
                    'segment mean':[np.mean(signal[simBorders[sim_i]:simBorders[sim_i+1]])],
                    'segment stdev':[np.std(signal[simBorders[sim_i]:simBorders[sim_i+1]])],
                    'segment length':[simBorders[sim_i+1] - simBorders[sim_i]]
                })

                segment_statistics = pd.concat([segment_statistics, base_sim_entry], ignore_index=True)

            for nano_i in np.random.choice(len(nanoBorders) - 1, size = 100, replace = False):#range(len(nanoBorders) - 1):

                nano_sim_entry = pd.DataFrame({
                    'segmentation type':['nanopolish'],
                    'read':[readid],
                    'segment mean':[np.mean(signal[nanoBorders[nano_i]:nanoBorders[nano_i+1]])],
                    'segment stdev':[np.std(signal[nanoBorders[nano_i]:nanoBorders[nano_i+1]])],
                    'segment length':[nanoBorders[nano_i+1] - nanoBorders[nano_i]]
                })

                segment_statistics = pd.concat([segment_statistics, nano_sim_entry], ignore_index=True)

    print('Done')
    print(f'Found {noSeg} reads without nanopolish segmentation')

    return read_statistics, segment_statistics

def plot(read_df : pd.DataFrame, segment_df : pd.DataFrame, path : str) -> None:
    '''
    Plot segmentation quality
    '''
    #  Read plots

    sns.set(style='whitegrid')
    sns.scatterplot(x='mean error', y='read segments stdev', data=read_df.loc[read_df['segmentation type'] == 'nanopolish'])
    ax = plt.gca()
    ax.set_title('Read mean error vs segment stdev')
    plt.tight_layout()
    plt.savefig(os.path.join(path, 'read_nano_std_error.png'))
    plt.close()
    
    g = sns.jointplot(data = read_df, x = 'read segments mean', y = 'read segments stdev', hue = 'segmentation type')
    g.fig.suptitle('Mean read segments means vs mean read segments stdevs\nin pA')
    g.fig.tight_layout()
    plt.savefig(os.path.join(path, 'read_segment_distributions.png'))
    plt.close()
    
    g = sns.jointplot(data = read_df.loc[read_df['segmentation type'] == 'nanopolish'], x = 'mean error', y = 'stdev error')
    g.fig.suptitle('Nanopolish segment mean error vs stdev error')
    g.fig.tight_layout()
    plt.savefig(os.path.join(path, 'read_nano_error.png'))
    plt.close()
    
    sns.histplot(data = read_df, x = 'basecall segmentation ratio', hue = 'segmentation type')
    ax = plt.gca()
    ax.set_title('Basecall-Segmentation ratio per segmentation type\n#-bases / #-segments')
    plt.tight_layout()
    plt.savefig(os.path.join(path, 'read_segmentation.png'))
    plt.close()
    
    # Segment plots
    g = sns.jointplot(data = segment_df, x = 'segment mean', y = 'segment stdev', hue = 'segmentation type')
    g.fig.suptitle('Segment mean vs segment stdev for all segments\n100 segments per read')
    g.fig.tight_layout()
    plt.savefig(os.path.join(path, 'segment_distributions.png'))
    plt.close()

    sns.histplot(data = segment_df, x = 'segment length', hue = 'segmentation type', bins = 60, stat = 'density')
    plt.title('Segment length distribution')
    plt.savefig(os.path.join(path, 'segment_length.png'))
    plt.close()

# TODO maybe change this to some kind of mapping, use bases/basecalling as index for segments
def getSegmentErrors(nanoBorders : np.ndarray, simBorders : np.ndarray) -> list:
    '''
    Return mean segmentation error

    Parameters
    ----------
    nanoBorders : sorted np.ndarray
    simBorders : sorted np.ndarray
    '''

    errors = []
    simIdx = 0

    for border in nanoBorders:
        while(border >= simBorders[simIdx] and simIdx < len(simBorders)):
            simIdx += 1
            
        ld = border - simBorders[simIdx - 1]
        rd = abs(simBorders[simIdx] - border)
        if ld <= rd:
            errors.append(ld)
        else:
            errors.append(rd)

    return errors

def readFastq(fq : str) -> dict:
    return SeqIO.to_dict(SeqIO.parse(open(fq),'fastq'))

def main() -> None:
    args = parse()
    sim5 = h5py.File(args.sim_fast5, 'r')
    nano5 = h5py.File(args.nanopolish_segmentation_fast5, 'r')
    fq_dict = readFastq(args.fastq)
    outpath = args.outpath

    if not os.path.exists(outpath):
        os.makedirs(outpath)

    print('Start comparison')
    read_df, segment_df = compare(sim5, nano5, fq_dict)
    plot(read_df, segment_df, outpath)

if __name__ == '__main__':
    main()
