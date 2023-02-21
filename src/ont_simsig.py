#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import os
import pandas as pd
from Bio import SeqIO
from numpy import Inf

from Simulator import RNASimulator
from Writer import RNAWriter

def parse() -> Namespace:

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter) 
    
    parser.add_argument('num_of_reads', type = int, help='Number of reads to simulate')
    parser.add_argument('outdir', type = str, help='Output directory to store files')
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-r', '--reference', default = None, type = str, help = 'Reference fasta file, currently only supporting a single reference')
    group.add_argument('-rl', '--reflen', default = None, type = int, help = 'Randomly generated reference length')

    parser.add_argument('--fullRef', default = False, action = 'store_true', help = 'Always draw full reference reads')
    parser.add_argument('--min_len', default = None, type = int, help = 'Minimum length of drawn reads')
    parser.add_argument('--max_len', default = None, type = int, help = 'Maximum length of drawn reads')
    parser.add_argument('--suffix', default = 'A'*50, help = 'Suffix sequence will be added to reference 3\' end')
    parser.add_argument('--model', type = str, default = os.path.join(os.path.dirname(__file__), '..', 'data', 'template_median69pA.model'), help = 'path to a .model file')
    parser.add_argument('--stdev_scale', type = float, default = 1.0, help = 'Scale to use for each models standard deviation.')
    parser.add_argument('--seed', default = None, type = int, help = 'Seed for the random generator')
    parser.add_argument('--segment_length', default = None, type = int, help = 'Set segment length for each kmer')
    parser.add_argument('--min_segment_length', default = 5, type = int, help = 'Minimum segment length')
    parser.add_argument('--max_segment_length', default = Inf, type = int, help = 'Maximum segment length')

    return parser.parse_args()

def readModel(model_path : str) -> dict:
    '''
    Reads the model csv

    Parameters
    ----------
    model_path : str
        path to model csv containing kmer, level_mean and level_stdv column
    
    Returns
    -------
    dict
        containing kmers as key and (mean, stdev) as tuple
    '''
    df = pd.read_csv(model_path, sep='\t')
    # kmers are stored from 3' to 5'
    return {key : (mean, stdev) for key, mean, stdev in zip(df['kmer'], df['level_mean'], df['level_stdv'])}

def readFasta(fa : str) -> tuple:
    fastas = SeqIO.to_dict(SeqIO.parse(open(fa),'fasta'))
    header = list(fastas.keys())[0]
    reference = list(fastas.keys())[1]
    return header, reference

def buildSimulator(model : dict, reference : str, reflen : int, suffix : str, seed : int, stdevScale : float, segmentLength : int, minL : int, maxL : float) -> tuple:

    if reference is not None:
        header, reference = readFasta(reference)
        rna = RNASimulator(model, reference=reference, suffix=suffix, seed=seed, stdevScale=stdevScale, setSegmentLength=segmentLength, shiftL = minL, maxL = maxL)
    else:
        rna = RNASimulator(model, refLength=reflen, suffix=suffix, seed=seed, stdevScale=stdevScale, setSegmentLength=segmentLength, shiftL = minL, maxL = maxL)
        reference = rna.getReference()
        header = None

    return header, reference, rna

def simulateReads(rnasimulator : RNASimulator, num_of_reads : int, fullRef : bool,  min_len : int, max_len : int) -> list:
    '''
    Returns
    signals : np.ndarray
        as [(signal, segments), ...]
        
        signal : np.ndarray
            a numpy array representing the simulated signal according the given kmer_model
        borders : np.ndarray
            an array containing the segment borders starting with 0
    '''

    if fullRef:
        return rnasimulator.drawRefSignals(num_of_reads)
    else:
        return rnasimulator.drawReadSignals(num_of_reads, minLen=min_len, maxLen=max_len)

def writeSignals(path : str, header : str, reference : str, signals) -> None:
    
    if not os.path.exists(path):
        os.makedirs(path)

    writer = RNAWriter(reference, path=path, tag = 'simulation', header = header)
    writer.writeReads(signals)

def main() -> None:
    args = parse()

    # required
    num_of_reads : int = args.num_of_reads
    outdir : str = args.outdir
    reference : str = args.reference
    reflen : int = args.reflen

    # optional
    min_len : int = args.min_len
    max_len : int = args.max_len
    suffix : str = args.suffix
    model : str = args.model
    fullRef : bool = args.fullRef
    stdev_scale : float = args.stdev_scale
    seed : int = args.seed
    segment_length : int = args.segment_length
    minL : int = args.min_segment_length
    maxL : float = args.max_segment_length

    print('Building reference ...')
    header, reference, rnasimulator = buildSimulator(readModel(model), reference, reflen, suffix, seed, stdev_scale, segment_length, minL, maxL)

    print('Simulating reads ...')
    signals = simulateReads(rnasimulator, num_of_reads, fullRef, min_len, max_len)

    print('Writing data ...')
    writeSignals(outdir, header, reference, signals)

    print('Done')

if __name__ == '__main__':
    main()