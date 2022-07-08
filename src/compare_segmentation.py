import h5py
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from Bio import SeqIO
import numpy as np
import pandas as pd
import seaborn as sns

def parse() -> Namespace:

    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
        ) 
    
    parser.add_argument('sim_fast5', type = str, help='FAST5 file with simulated reads')
    parser.add_argument('nanopolish_segmentation_fast5', type = str, help='FAST5 file with nanopolish segmentation')
    parser.add_argument('fastq', type = str, help = 'Basecalled reads')

    return parser.parse_args()

def compare(sim5 : h5py.File, nano5 : h5py.File, fastqs : dict):
    
    # store key features for each read describing the segmentation and basecalling quality
    # TODO change to pandas dataframe for plotting
    segmentation_quality = {}

    for readid in nano5:

        # key in sim5, ONT FAST5 format
        read = f'read_{readid}'

        nanoBorders = nano5[readid][:, 0]
        simBorders = sim5[read]['Raw/Borders'][:]
        nSimSegments = sim5[read]['Raw'].attrs['num_segments'] # ~ #bases: #bases = #segments + 4, #segments = #5mers
        signal = sim5[read]['Raw/Signal'][:]
        assert (nSimSegments + 1) == len(simBorders)

        read : str = str(fastqs[readid].seq)
        # phred : list[int] = fastqs[readid].letter_annotations['phred_quality']

        segmentation_quality[readid]['target_mean'] = np.mean([np.mean(signal[simBorders[i]:simBorders[i+1]]) for i in range(nSimSegments)])
        # intrasegmental standard deviation
        segmentation_quality[readid]['target_stdev'] = np.mean([np.std(signal[simBorders[i]:simBorders[i+1]]) for i in range(nSimSegments)])

        segmentation_quality[readid]['nano_mean'] = np.mean([np.mean(signal[nanoBorders[i]:nanoBorders[i+1]]) for i in range(len(nanoBorders) - 1)])
        segmentation_quality[readid]['nano_stdev'] = np.mean([np.std(signal[nanoBorders[i]:nanoBorders[i+1]]) for i in range(len(nanoBorders) - 1)])

        segmentation_quality[readid]['target_fraction'] = len(read)/nSimSegments
        segmentation_quality[readid]['nano_fraction'] = len(read)/(len(simBorders) - 1)

        segmentation_quality[readid]['nano_error'] = getSegmentError(nanoBorders, simBorders)

def getSegmentError(nanoBorders : np.ndarray, simBorders : np.ndarray) -> float:
    '''
    Return mean segmentation error

    Parameters
    ----------
    nanoBorders : sorted np.ndarray
    simBorders : sorted np.ndarray
    '''
    error = []
    simIdx = 0

    for border in nanoBorders:
        while(border <= simBorders[simIdx]):
            simIdx += 1
        error.append(np.abs(border - simBorders[simIdx - 1]))

    return np.mean(error)

def readFastq(fq : str) -> dict:
    return SeqIO.to_dict(SeqIO.parse(open(fq),'fastq'))

def main() -> None:
    args = parse()
    sim5 = h5py.File(args.sim_fast5, 'r')
    nano5 = h5py.File(args.nanopolish_segmentation_fast5, 'r')
    fq_dict = readFastq(args.fastq)

    compare(sim5, nano5, fq_dict)

if __name__ == '__main__':
    main()