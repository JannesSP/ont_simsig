from typing import Iterable, Tuple
import numpy as np

class RNASimulator():

    def __init__(self, kmer_model : dict, length : int = 2000, reference : str = None, alphabet : str = 'ACGT', checkNucl : bool = True, prefix : str = '', suffix : str = '', seed : int = None, stdev_scale : float = 1.0, set_segment_length : int = None) -> None:
        '''
        Parameters
        ----------
        kmer_model : dict
            Keys are RNA 5mers (``3' to 5' orientation``) and values are tuples (mean, stdev) of the gaussian signal distribution
        length : int
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
        '''
        if checkNucl:
            for nucl in alphabet:
                assert nucl in 'ACGT', f'Nucleotide {nucl} not part of ACGT in alphabet'
            for nucl in prefix:
                assert nucl in 'ACGT', f'Nucleotide {nucl} not part of ACGT in prefix'
            for nucl in suffix:
                assert nucl in 'ACGT', f'Nucleotide {nucl} not part of ACGT in suffix'

        assert 0.0 <= stdev_scale <= 1.0, 'Standard deviation scale must be in [0.0, 1.0]'
        if set_segment_length is not None:
            assert set_segment_length >= 5

        self.alphabet = alphabet
        self.stdev_scale = stdev_scale
        self.set_segment_length = set_segment_length
        self.npAlphabet = np.array(list(alphabet))
        self.kmer_model = kmer_model
        self.__setSeed(seed)

        # mean for exponential distribution
        self.exp_m = 55.7

        # special sequence that should always appear at 5' end of sequence like adapters, barcode
        self.prefix = prefix
        # special sequence that should always appear at 3' end of sequence like adapters, barcode, polyA
        self.suffix = suffix
        self.__generate(reference, length)
        self.simulatedReads = 0

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

    def __generate(self, reference, length) -> None:
        '''
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
            self.length = len(self.prefix) + length + len(self.suffix)
        else:
            assert len(reference) > 4, f'Reference sequence too small ({len(reference)}), cannot initialize'
            self.reference = self.prefix + reference + self.suffix
            self.length = len(self.reference)

    def __drawSegmentLength(self) -> int:
        '''
        Model to simulate the segment length in a RNA read
        '''
        if self.set_segment_length is not None:
            return self.set_segment_length
        else:
            return np.random.exponential(self.exp_m, 1).astype(int).item() + 5 # minimum segment length I saw in nanopolish segmentation was 5

    def __drawSegmentSignal(self, kmer : str, length : int) -> np.ndarray:
        assert len(kmer) == 5, f'Kmer must have length 5 not {len(kmer)}'
        segment=np.random.normal(self.kmer_model[kmer][0], self.kmer_model[kmer][1]*self.stdev_scale, length)
        # print(self.kmer_model[kmer][0], self.kmer_model[kmer][1]*self.stdev_scale, segment)
        return segment

    def getNumSimReads(self) -> int:
        return self.simulatedReads

    def getRefLength(self) -> int:
        return self.length

    def getReference(self) -> str:
        '''
        Returns
        -------
        reference : str
            reference string in 5' to 3' orientation
        '''
        return self.reference

    def drawReadSignals(self, n : int, max_len : int, min_len : int = 5) -> Iterable[Tuple[np.ndarray, np.ndarray]]:
        '''
        Generates n read signal starting from 3' (RNA) end of the reference and stopping after a uniformly drawn number of bases of the interval [min_len, max_len)

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
        sim_signals : np.ndarray
            a numpy array representing the simulated signal according the given kmer_model
        borders : np.ndarray
            an array containing the segment borders starting with 0
        '''
        assert n > 0
        assert min_len < max_len
        assert max_len < self.length
        assert min_len >= 5
        self.simulatedReads += n
        return np.array([self.__drawSignal(stop = np.random.randint(min_len, max_len, size = 1, dtype = int).item()) for _ in range(n)])

    def drawReadSignal(self, max_len : int, min_len : int = 5) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Generates n read signal starting from 3' (RNA) end of the reference and stopping after a uniformly drawn number of bases

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
        return self.__drawSignal(stop = np.random.randint(min_len, max_len, size = 1, dtype = int).item())

    def drawRefSignals(self, n : int) -> Iterable[Tuple[np.ndarray, np.ndarray]]:
        '''
        Generates n read signal starting from 3' (RNA) end of the reference and stopping after a uniformly drawn number of bases

        Parameters
        ----------
        n : int
            number of generated signals
        '''
        assert n > 0
        sims = []
        for i in range(n):
            if (i+1)%10==0:
                print(f'Simulating read {i + 1}\\{n} ...', end = '\r')
            sims.append(self.__drawSignal())
        print(f'Done simulating {n} reads  ')
        self.simulatedReads += n
        return sims

    def drawRefSignal(self) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Generates 1 read signal starting from 3' (RNA) end of the reference and stopping after a uniformly drawn number of bases
        '''
        self.simulatedReads += 1
        return self.__drawSignal()

    def __drawSignal(self, stop : int = None) -> Tuple[np.ndarray, np.ndarray]:
        '''
        Generates a read signal from the whole reference

        Parameters
        ----------
        stop : int
            simulate signal from the 3' end to this base position (base included and 1-based)

        Returns
        -------
        sim_signals : np.ndarray
            a numpy array representing the simulated signal according the given kmer_model from 5' to 3' end
        borders : np.ndarray
            an array containing the segment borders starting with 0 at the 5' end
        '''
        
        # change orientation of the reference to 3' to 5' (RNA is sequenced from 3' end)
        if stop is not None:
            reference = self.reference[::-1][:stop]
        else:
            reference = self.reference[::-1]

        assert len(reference) > 4, f'Reference sequence too small ({len(reference)}) for simulation'

        # current length of the simulated signal
        signal_pointer = 0
        init_len = len(reference) * int(self.exp_m)
        sim_signal = np.zeros(init_len, dtype = float)

        border_pionter = 1
        borders = np.zeros(len(reference) - 3, dtype = int)

        for n in range(len(reference) - 4):
            # kmer_model is in 3' -> 5' orientation, same as reference here
            kmer = reference[n:n+5]
            segment_length = self.__drawSegmentLength()

            if signal_pointer + segment_length >= len(sim_signal):
                sim_signal = np.append(sim_signal, np.zeros(init_len, dtype = float))

            sim_signal[signal_pointer : signal_pointer + segment_length] = self.__drawSegmentSignal(kmer, segment_length)
            signal_pointer += segment_length

            borders[border_pionter] = signal_pointer
            border_pionter += 1

        return sim_signal[:signal_pointer], borders[:border_pionter]