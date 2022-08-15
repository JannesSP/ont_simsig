# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from typing import List, Tuple
import pysam
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from Bio import SeqIO
import matplotlib.pyplot as plt
import os
import numpy as np
import seaborn as sns

def parse() -> Namespace:

    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
        ) 
    
    parser.add_argument('bam_file', type = str, help='bam file input')
    parser.add_argument('reference_fasta', type = str, help='Fasta reference file')

    return parser.parse_args()

def readFasta(fasta : str) -> Tuple[str, str]:
    '''
    Reads only first fasta sequence in the .FASTA file

    Parameters
    ----------
    fasta : str
        Path to fasta file
    '''
    fasta_sequences = SeqIO.parse(open(fasta),'fasta').__next__()
    print(f'Reading fasta:\n{fasta_sequences.id}\n{str(fasta_sequences.seq[:100])}...')
    return fasta_sequences.id, str(fasta_sequences.seq)

def countMapping(bam : str, ref_header : str, ref_sequence : str) -> List[dict]:

    # counts of bases per ref position, N == total number of bases mapped at this position
    bases = []
    idx = 0

    samfile = pysam.AlignmentFile(bam, 'rb')
    for idx, pileupcolumn in enumerate(samfile.pileup(ref_header, 1, len(ref_sequence))):

        if (idx + 1)%100 == 0:
            print(f'Getting data {idx + 1}\{len(ref_sequence)}', end = '\r')
        
        bases.append(
            {'A':0,
            'C':0,
            'G':0,
            'T':0,
            'N':pileupcolumn.nsegments,
            'counts':0,
            'pos':pileupcolumn.reference_pos,
            'ref_base':ref_sequence[pileupcolumn.reference_pos],
            'del':0,
            'skip':0
            })

        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del:
                bases[idx]['del'] += 1
            elif pileupread.is_refskip:
                bases[idx]['skip'] += 1
            else:
                # query position is None if is_del or is_refskip is set.
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                # print(pileupread.alignment.query_sequence[pileupread.query_position], ref_sequence[pileupcolumn.reference_pos], pileupread.alignment.query_sequence[pileupread.query_position] == ref_sequence[pileupcolumn.reference_pos])
                bases[idx][base] += 1
                bases[idx]['counts'] += 1

    print(f'Found alignments for {idx} positions')

    # print(f'Got data for {idx}\{len(ref_sequence)} positions')
    return bases

def plotMapping(x : np.ndarray, id : np.ndarray, cov : np.ndarray, d : np.ndarray, path : str = None) -> None:

    sum = open(os.path.join(path, 'counts.csv'), 'w+')
    sum.write('feature,mean,stdev,median,min,max\n')
    sum.write(f'identity,{np.mean(id)},{np.median(id)},{np.min(id)},{np.max(id)}\n')
    sum.write(f'coverage,{np.mean(cov)},{np.median(cov)},{np.min(cov)},{np.max(cov)}\n')
    sum.write(f'deletions,{np.mean(d)},{np.median(d)},{np.min(d)},{np.max(d)}\n')

    assert len(x) == len(id)

    plt.figure(figsize = (8,6), dpi=800)
    lns1=plt.plot(x, id, label = f'identity rate mean={np.mean(id):.3f}', alpha = 0.5, color='blue')
    lns2=[plt.hlines(np.mean(id), min(x), max(x), linestyles='--', color='blue', label=f'mean={np.mean(id):.3f}')]
    lns3=plt.plot(x, d, label = f'deletion rate mean={np.mean(d):.2f}', color = 'red', alpha = 0.5)
    lns4=[plt.hlines(np.mean(d), min(x), max(x), linestyles='--', color='red', label=f'mean={np.mean(d):.2f}')]
    plt.title('basecalls result')
    plt.ylabel('identity and deletion rate')
    plt.xlabel('reference position')
    plt.grid(True, 'both', 'both', alpha = 0.8)
    ax2 = plt.twinx()
    lns5=ax2.plot(x, cov, label = f'coverage mean={np.mean(cov):.2f}', color = 'orange')
    ax2.set_ylabel('coverage')
    lns = lns1+lns2+lns3+lns4+lns5
    labs = [l.get_label() for l in lns]
    plt.legend(lns, labs, loc=0)
    plt.tight_layout()

    if path is None:
        plt.show()
    else:
        plt.savefig(os.path.join(path, 'basecall_results.png'))

    plt.close()

    # sort for Identity
    sortIds = np.sort(id)
    sortDels = np.sort(d)

    plt.figure(figsize = (8,6), dpi=800)
    lns1=plt.plot(sortIds, label = f'sorted identity rate', color='blue')
    lns2=[plt.hlines(np.mean(id), min(x), max(x), linestyles='--', color='blue', label=f'mean={np.mean(id):.3f}')]
    plt.title('basecalls result sorted for identity')
    plt.ylabel('identity')
    plt.tick_params(axis='y', which='both', labelright=True)
    plt.grid(True, 'both', 'both')
    lns3=plt.plot(sortDels, label = f'sorted deletion rate', color = 'red')
    lns4=[plt.hlines(np.mean(d), min(x), max(x), linestyles='--', color='red', label=f'mean={np.mean(d):.2f}')]
    lns = lns1+lns2+lns3+lns4
    labs = [l.get_label() for l in lns]
    plt.legend(lns, labs, loc=0)
    plt.tight_layout()

    if path is None:
        plt.show()
    else:
        plt.savefig(os.path.join(path, 'basecall_results_sorted.png'))

    plt.close()

    # identity histplot

    sns.histplot(x = sortIds, kde=True, stat="density", bins = 100)
    plt.xlabel('Identity')
    plt.ylabel('Density')
    plt.title('Identity distribution of simulated data')
    plt.tight_layout()

    if path is None:
        plt.show()
    else:
        plt.savefig(os.path.join(path, 'basecall_identity.png'))

    plt.close()

def main() -> None:
    args = parse()
    bam = args.bam_file
    fasta = args.reference_fasta

    savedir = os.path.dirname(bam)

    counts = countMapping(bam, *readFasta(fasta))

    id = np.array([(pos[pos['ref_base']]/pos['counts']) if pos['counts'] != 0 else 0 for pos in counts])
    cov = np.array([pos['N'] for pos in counts])
    x = np.array([pos['pos'] for pos in counts])
    d = np.array([pos['del']/pos['N'] for pos in counts])

    assert len(x) > 0, f'No mapped reads found to plot'
    
    plotMapping(x, id, cov, d, path = savedir)

if __name__ == '__main__':
    main()

