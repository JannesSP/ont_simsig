# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import os
from Bio import SeqIO
from numpy import Inf

from __init__ import __version__
from Simulator import RNASimulator
from Writer import RNAWriter

def parse() -> Namespace:

    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter) 
    
    parser.add_argument('numOfReads', type = int, help='Number of reads to simulate')
    parser.add_argument('outdir', type = str, help='Output directory to store files')
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-r', '--reference', default = None, type = str, help = 'Reference fasta file, currently only supporting a single reference')
    group.add_argument('-rl', '--refLen', default = None, type = int, help = 'Randomly generated reference of given length')
    parser.add_argument('--fullRef', default = False, action = 'store_true', help = 'Always draw full reference reads')
    parser.add_argument('--minReadLen', default = 100, type = int, help = 'Minimum length of drawn reads')
    parser.add_argument('--maxReadLen', default = None, type = int, help = 'Maximum length of drawn reads')
    parser.add_argument('--suffix', default = 'A'*50, help = 'Suffix sequence will be added to reference 3\' end')
    parser.add_argument('--stdevScale', type = float, default = 1.0, help = 'Scale to use for each models standard deviation.')
    parser.add_argument('--seed', default = None, type = int, help = 'Seed for the random generator')
    parser.add_argument('--segLen', default = None, type = int, help = 'Set segment length for each kmer')
    parser.add_argument('--maxSegLen', default = Inf, type = int, help = 'Maximum segment length')

    return parser.parse_args()

def readFasta(fa : str) -> tuple:
    fastas = SeqIO.to_dict(SeqIO.parse(open(fa),'fasta'))
    header = list(fastas.keys())[0]
    return header, fastas[header].seq

def buildSimulator(reference : str, reflen : int, suffix : str, seed : int, stdevScale : float, segmentLength : int, maxL : float) -> tuple:

    if reference is not None:
        header, reference = readFasta(reference)
        rna = RNASimulator(reference=reference, suffix=suffix, seed=seed, stdevScale=stdevScale, setSegmentLength=segmentLength, maxL = maxL)
    else:
        rna = RNASimulator(refLength=reflen, suffix=suffix, seed=seed, stdevScale=stdevScale, setSegmentLength=segmentLength, maxL = maxL)
        reference = rna.getReference()
        header = None

    return header, reference, rna

def simulateReads(rnasimulator : RNASimulator, numOfReads : int, fullRef : bool,  minReadLen : int, maxReadLen : int) -> list:
    '''
    Returns
    signals : np.ndarray
        as [(signal, segments, read), ...]
        
        signal : np.ndarray
            a numpy array representing the simulated signal according the given kmer_model
        borders : np.ndarray
            an array containing the segment borders starting with 0
    '''
    if fullRef:
        return rnasimulator.drawRefSignals(numOfReads)
    else:
        return rnasimulator.drawReadSignals(numOfReads, minLen=minReadLen, maxLen=maxReadLen)

def writeSignals(path : str, header : str, reference : str, signals) -> None:
    
    if not os.path.exists(path):
        os.makedirs(path)

    writer = RNAWriter(reference, path=path, tag='simulation', header=header)
    writer.writeReads(signals)

def main() -> None:
    args = parse()

    # required
    numOfReads : int = args.numOfReads
    outdir : str = args.outdir
    reference : str = args.reference
    refLen : int = args.refLen

    # optional
    minReadLen : int = args.minReadLen
    maxReadLen : int = args.maxReadLen
    if maxReadLen is None:
        maxReadLen = refLen
    suffix : str = args.suffix
    fullRef : bool = args.fullRef
    stdevScale : float = args.stdevScale
    seed : int = args.seed
    segmentLength : int = args.segLen
    maxL : float = args.maxSegLen

    print('Building reference ...')
    header, reference, rnasimulator = buildSimulator(reference, refLen, suffix, seed, stdevScale, segmentLength, maxL)
    print('Simulating reads ...')
    results = simulateReads(rnasimulator, numOfReads, fullRef, minReadLen, maxReadLen)
    print('Writing data ...')
    writeSignals(outdir, header, reference, results)
    print('Done')

if __name__ == '__main__':
    main()