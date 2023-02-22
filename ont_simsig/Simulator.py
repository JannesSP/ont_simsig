# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from typing import Iterable, Tuple
import numpy as np
from os.path import join

class RNASimulator():

    def __init__(
        self,
        refLength : int = 2000,
        reference : str = None,
        alphabet : str = 'ACGT',
        checkNucl : bool = True,
        prefix : str = '',
        suffix : str = '',
        seed : int = None,
        stdevScale : float = 1.0,
        setSegmentLength : int = None,
        # exp_m : float = 31.2, old length simulation of was with exponential distribution
        # shiftL : int = 5, old shift is now contained in the kmer length distribution file
        maxL : int = np.inf) -> None:

        '''
        Parameters
        ----------
        ref_length : int
            Length for the randomly generated reference
        reference : str
            Explicitly set reference for simulation
        alphabet : str
            Set of nucleotides used for reference generation
        checkNucl : bool
            Check alphabet, prefix and suffix for ACGT
        prefix : str
            Nucleotide sequence added to 5' end
        suffix : str
            Nucleotide sequence added to 3' end
        seed : int
            Set seed for the random generators
        maxL : int
            Maximum segment length, default np.inf
        '''
        if checkNucl:
            for nucl in alphabet:
                assert nucl in 'ACGT', f'Nucleotide {nucl} not part of ACGT in alphabet'
            for nucl in prefix:
                assert nucl in 'ACGT', f'Nucleotide {nucl} not part of ACGT in prefix'
            for nucl in suffix:
                assert nucl in 'ACGT', f'Nucleotide {nucl} not part of ACGT in suffix'

        assert 0.0 <= stdevScale <= 1.0, 'Standard deviation scale must be in [0.0, 1.0]'
        if setSegmentLength is not None:
            assert setSegmentLength >= 5

        self.alphabet = alphabet
        self.stdevScale = stdevScale
        self.setSegmentLength = setSegmentLength
        self.npAlphabet = np.array(list(alphabet))
        self.__setSeed(seed)
        self.maxL = maxL
        self.prefix = prefix # special sequence that should always appear at 5' end of sequence like adapters, barcode
        self.suffix = suffix # special sequence that should always appear at 3' end of sequence like adapters, barcode, polyA
        self.__generateReference(reference, refLength)
        self.simulatedReads = 0

        # load distributions models
        self.__loadKmerModels()
        self.__loadLenModels()

    def __loadKmerModels(self) -> None:
        '''
        Loads signal distirbution kmer models from template_median69pA.model
        '''
        self.sigModels = {}
        with open(join(__file__, '..', '..', 'data', 'template_median69pA.model'), 'r') as models:
            models.readline() # skip header
            for line in models:
                kmer, mean, stdev, _, _, _, _ = line.strip().split('\t')
                self.sigModels[kmer[::-1]] = (float(mean), float(stdev)) # kmers are in 3'->5' orientation in template file from ONT

    def __loadLenModels(self) -> None:
        '''
        Loads segment length kmer models from kmer_nbin.csv into a dictionary.
        '''
        self.lenModels = {}
        with open(join(__file__, '..', '..', 'data', 'kmer_nbin.csv'), 'r') as models:
            models.readline() # skip header
            for line in models:
                _, kmer, p, r, shift, _ = line.strip().split(',')
                self.lenModels[kmer] = (float(p), float(r), int(shift)) # kmers are in 5'->3' orientation in our segment length file

    def __setSeed(self, seed : int) -> None:
        '''
        Set the numpy seed

        Parameters
        ----------
        seed : int
            Seed to set with np.random.seed(seed)
        '''
        if seed is not None:
            np.random.seed(seed)
            self.seed = seed

    def __generateReference(self, reference, length) -> None:
        '''
        Generate the reference sequence in 5'->3' orientation

        Parameters
        ----------
        reference : str
            Reference string to use for simulation
        length : int
            Length of randomly generated reference for simulation
        '''

        print('Generating reference')

        if reference is None:
            self.reference = self.prefix + ''.join(np.random.choice(self.npAlphabet, size = length)) + self.suffix
            self.refLength = len(self.prefix) + length + len(self.suffix)
        else:
            assert len(reference) > 4, f'Reference sequence too small ({len(reference)}), cannot initialize'
            self.reference = self.prefix + reference + self.suffix
            self.refLength = len(self.reference)

    def __drawSegmentLength(self, kmer : str) -> int:
        '''
        Model to simulate the segment length in a RNA read.
        Will return set segment length, if set.
        Otherwise returns random integer segment length from a negative binomial distribution according to kmer_nbin.csv.
        
        Parameters
        ----------
        kmer : str
            5 mer RNA sequence in 5'->3' orientation

        Returns
        -------
        segment_length : int
            Segment length for the given kmer
        '''

        if self.setSegmentLength is not None:
            return self.setSegmentLength
        else:
            # return min(np.random.exponential(self.exp_m, 1).astype(int).item() + self.shiftL, self.maxL) # minimum segment length I saw in nanopolish segmentation was 5
            r, n, shiftL = self.lenModels[kmer]
            return min(np.random.negative_binomial(n, r, 1).item() + shiftL, self.maxL)

    def __drawSegmentSignal(self, kmer : str, length : int) -> np.ndarray:
        mean, stdev = self.sigModels[kmer]
        return np.random.normal(mean, stdev*self.stdevScale, length)

    def getNumSimReads(self) -> int:
        return self.simulatedReads

    def getRefLength(self) -> int:
        return self.refLength

    def getReference(self) -> str:
        '''
        Returns
        -------
        reference : str
            reference string in 5' to 3' orientation
        '''
        return self.reference

    def drawReadSignals(self, n : int, maxLen : int, minLen : int = 5) -> Iterable[Tuple[np.ndarray, np.ndarray]]:
        '''
        Generates n read signals starting from 3' (RNA) end of the reference and stopping after a uniformly drawn number of bases of the interval [min_len, max_len)

        Parameters
        ----------
        n : int
            number of generated signals
        max_len : int
            max length of the used reference
        min_len : int
            min length of the used reference

        Returns
        -------
        signals : np.ndarray
            as [(signal, segments), ...]
            
            signal : np.ndarray
                a numpy array representing the simulated signal according the given kmer_model
            borders : np.ndarray
                an array containing the segment borders starting with 0
        '''
        assert n > 0
        assert minLen is not None
        assert minLen < maxLen
        assert maxLen < self.refLength
        assert minLen >= 5
        self.simulatedReads += n
        return np.array([self.__drawSignal(stop = np.random.randint(minLen, maxLen, size = 1, dtype = int).item()) for _ in range(n)])

    def drawReadSignal(self, maxLen : int, minLen : int = 5) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Generates 1 read signal starting from 3' (RNA) end of the reference and stopping after a uniformly drawn number of bases

        Parameters
        ----------
        max_len : int
            max length of the used reference
        min_len : int
            min length of the used reference

        Returns
        -------
        sim_signals : np.ndarray
            a numpy array representing the simulated signal according the given kmer_model
        borders : np.ndarray
            an array containing the segment borders starting with 0
        '''
        self.simulatedReads += 1
        return self.__drawSignal(stop = np.random.randint(minLen, maxLen, size = 1, dtype = int).item())

    def drawRefSignals(self, n : int, segmentLengths : Iterable[Iterable[int]] = None) -> Iterable[Tuple[np.ndarray, np.ndarray]]:
        '''
        Generates n read signal starting from 3' (RNA) end of the reference and stopping after a uniformly drawn number of bases

        Parameters
        ----------
        n : int
            number of generated read signals
        segment_lengths : Iterable[Iterable[int]]
            Array of integer arrays containing the number of signal data points for each 5mer segment (= #bases - 4) for each simulated read (n)

        Returns
        -------
        signals : np.ndarray
            as [(signal, segments), ...]
            
            signal : np.ndarray
                a numpy array representing the simulated signal according the given kmer_model
            borders : np.ndarray
                an array containing the segment borders starting with 0
        '''
        assert n > 0
        sims = []
        for i in range(n):
            if (i+1)%10==0:
                print(f'Simulating read {i + 1}\\{n} ...', end = '\r')
            if segmentLengths is not None:
                sims.append(self.__drawSignal(segmentLengths=segmentLengths[i]))
            else:
                sims.append(self.__drawSignal())
        print(f'Done simulating {n} reads  ')
        self.simulatedReads += n
        return sims

    def drawRefSignal(self, segmentLengths : Iterable[int] = None) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Generates 1 read signal starting from 3' (RNA) end of the reference and stopping after a uniformly drawn number of bases

        Parameters
        ----------
        segment_lengths : Iterable[int]
            Integer array with the number of signal data points for each 5mer segment (= #bases - 4)
        '''
        self.simulatedReads += 1
        return self.__drawSignal(segmentLengths = segmentLengths)

    def __drawTransition(x : float, y : float) -> np.array:
        '''
        Generates transition datapoints between segments
        '''
        pass

    def __drawSignal(self, stop : int = None, segmentLengths : Iterable[int] = None) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Generates a read signal from the whole reference

        Parameters
        ----------
        stop : int
            simulate signal from the 3' end to this base position (base included and 1-based)
        segment_lengths : Iterable[int]
            Integer array with the number of signal data points for each 5mer segment (= #bases - 4)

        Returns
        -------
        sim_signals : np.ndarray
            a numpy array representing the simulated signal in 3'->5' orientation
        borders : np.ndarray
            an array containing the segment borders starting with 0 from 3' end
        '''
        
        # change orientation of the reference to 3' to 5' (RNA is sequenced from 3' end)
        if stop is not None:
            reference = self.reference[::-1][:stop]
        else:
            reference = self.reference[::-1]
            stop = len(reference)

        if segmentLengths is not None:
            assert stop == len(segmentLengths) - 4, 'Len of segment lenghts list must be same as number of 5mer signals to draw'

        assert len(reference) > 4, f'Reference sequence too small ({len(reference)}) for simulation'

        # current length of the simulated signal
        signalPointer = 0
        initLen = len(reference) * int(35) # 35.163636... is the mean of the negative binomial for 'all' kmers
        simSignal = np.zeros(initLen, dtype = float)

        borderPionter = 1
        borders = np.zeros(len(reference) - 3, dtype = int)

        # loop over reference in 3'->5' orientation
        for i in range(len(reference) - 4):
            kmer = reference[i:i+5][::-1] # model kmers are in 5'->3' direction, change it to 3'->5' for simulation of RNA signal
            segmentLength = segmentLengths[i] if segmentLengths is not None else self.__drawSegmentLength(kmer)

            if signalPointer + segmentLength >= len(simSignal):
                simSignal = np.append(simSignal, np.zeros(initLen, dtype = float))

            simSignal[signalPointer : signalPointer + segmentLength] = self.__drawSegmentSignal(kmer, segmentLength)
            signalPointer += segmentLength

            borders[borderPionter] = signalPointer
            borderPionter += 1

        return simSignal[:signalPointer], borders[:borderPionter]