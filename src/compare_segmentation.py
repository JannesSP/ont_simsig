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

    return parser.parse_args()

def compare(sim5 : h5py.File, nano5 : h5py.File, fastqs : dict) -> pd.DataFrame:
    
    # store key features for each read describing the segmentation and basecalling quality
    df = pd.DataFrame(
        columns=['type', 'segment_mean', 'segment_stdev', 'segmentation', 'mean_error', 'stdev_error', 'max_error', 'min_error'],
        dtype=['str', 'float','float','float','float','float','float','float'])

    for readid in nano5:

        # key in sim5, ONT FAST5 format
        read = f'read_{readid}'

        nanoBorders = nano5[readid][:, 0]
        simBorders = sim5[read]['Raw/Borders'][:]
        nSimSegments = sim5[read]['Raw'].attrs['num_segments'] # ~ #bases: #bases = #segments + 4, #segments = #5mers
        signal = sim5[read]['Raw/Signal'][:]
        assert (nSimSegments + 1) == len(simBorders)
        e = getSegmentErrors(nanoBorders, simBorders)
        read : str = str(fastqs[readid].seq)
        # phred : list[int] = fastqs[readid].letter_annotations['phred_quality']

        sim_entry = pd.DataFrame({
            'type':['simulation'],
            'segment_mean':[np.mean([np.mean(signal[simBorders[i]:simBorders[i+1]]) for i in range(nSimSegments)])],
            'segment_stdev':[np.mean([np.std(signal[simBorders[i]:simBorders[i+1]]) for i in range(nSimSegments)])],
            'segmentation':[len(read)/nSimSegments],
            'mean_error':[np.nan],
            'stdev_error':[np.nan],
            'max_error':[np.nan],
            'min_error':[np.nan]
        })

        nano_entry = pd.DataFrame({
            'type':['nanopolish'],
            'segment_mean':[np.mean([np.mean(signal[nanoBorders[i]:nanoBorders[i+1]]) for i in range(len(nanoBorders) - 1)])],
            'segment_stdev':[np.mean([np.std(signal[nanoBorders[i]:nanoBorders[i+1]]) for i in range(len(nanoBorders) - 1)])],
            'segmentation':[len(read)/(len(simBorders) - 1)],
            'mean_error':[np.mean(e)],
            'stdev_error':[np.std(e)],
            'max_error':[np.max(e)],
            'min_error':[np.min(e)]
        })
        
        df = df.append([df, sim_entry, nano_entry], ignore_index=True)
        
    return df

def plot(df : pd.DataFrame, path : str) -> None:
    '''
    Plot segmentation quality
    '''
    sns.set(style='whitegrid')
    sns.lmplot(x='segment_stdev', y='mean_error', data=df.loc[df['type'] == 'nanopolish'])
    plt.savefig(os.path.join(path, 'nano_std_error.png'))
    plt.close()
    
    sns.jointplot(data = df, x = 'segment_mean', y = 'segment_stdev', hue = 'type')
    plt.savefig(os.path.join(path, 'segment_distritbuions.png'))
    plt.close()
    
    sns.jointplot(data = df.loc[df['type'] == 'nanopolish'], x = 'mean_error', y = 'stdev_error')
    plt.savefig(os.path.join(path, 'nano_error.png'))
    plt.close()
    
    sns.histplot(data = df, x = 'segmentation', hue = 'type')
    plt.savefig(os.path.join(path, 'segmentation.png'))
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
        while(border >= simBorders[simIdx]):
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

    df = compare(sim5, nano5, fq_dict)

if __name__ == '__main__':
    main()
