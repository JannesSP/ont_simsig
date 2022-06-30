import random
import numpy as np

class RNASimulator():

    def __init__(self, kmer_model : dict, length : int = 2000, reference : str = None, alphabet : str = 'ACGT', prefix : str = '', suffix : str = '') -> None:
        '''
        Parameters
        ----------
        kmer_model : dict
            Keys should be 5mers and values should be a tuple containing mean and stdev of the signal distribution
        length : int
            Length for the randomly generated reference
        reference : str
            Explicitly set reference for simulation
        alphabet : str
            Set of nucleotides used for reference generation
        prefix : str
            Nucleotide sequence added to 5' end
        suffix : str
            Nucleotide sequence added to 3' end
        '''
        self.alphabet = alphabet
        self.kmer_model = kmer_model

        # mean for exponential distribution
        self.exp_m = 55.7

        # special sequence that should always appear at 5' end of sequence like adapters, barcode
        self.prefix = prefix
        # special sequence that should always appear at 3' end of sequence like adapters, barcode, polyA
        self.suffix = suffix
        self.__generate__(reference, length)

    def __generate__(self, reference, length) -> None:
        '''
        Parameters
        ----------
        reference : str
            Reference string to use for simulation
        length : int
            Length of randomly generated reference for simulation
        '''

        if reference is None:
            self.reference = self.prefix + ''.join(random.choices(self.alphabet, k=length)) + self.suffix
            self.length = len(self.prefix) + length + len(self.suffix)
        else:
            assert len(reference) > 4, f'Reference sequence too small ({len(reference)}), cannot initialize'
            self.reference = self.prefix + reference + self.suffix
            self.length = len(self.reference)

        for nucl in self.reference:
            assert nucl in self.alphabet, f'Reference contains a nucleotide ({nucl}) that is not in the alphabet ({self.alphabet})!'

    def __drawSegmentLength__(self) -> int:
        '''
        Model to simulate the segment length in a RNA read
        '''
        return np.random.exponential(self.exp_m, 1).astype(int).item()

    def __drawSegmentSignal__(self, kmer : str, length : int) -> np.ndarray:
        assert len(kmer) == 5, f'Kmer must have length 5 not {len(kmer)}'
        return np.random.normal(self.kmer_model[kmer][0], self.kmer_model[kmer][1], length)

    def getRefLength(self) -> int:
        return self.length

    def getReference(self) -> str:
        return self.reference

    def drawReadSignals(self, n : int, max_len : int, min_len : int = 5) -> list[tuple[np.ndarray, np.ndarray]]:
        '''
        Generates n read signal starting from 3' (RNA) end of the reference and stopping after a uniformly drawn number of bases

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
        return [self.__drawSignal__(stop = np.random.randint(min_len, max_len, size = 1, dtype = int).item()) for _ in range(n)]

    def drawReadSignal(self, max_len : int, min_len : int = 5) -> tuple[np.ndarray, np.ndarray]:
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
        return self.__drawSignal__(stop = np.random.randint(min_len, max_len, size = 1, dtype = int).item())

    def drawRefSignals(self, n : int) -> list[tuple[np.ndarray, np.ndarray]]:
        '''
        Generates n read signal starting from 3' (RNA) end of the reference and stopping after a uniformly drawn number of bases

        Parameters
        ----------
        n : int
            number of generated signals
        '''
        assert n > 0
        return [self.__drawSignal__() for _ in range(n)]

    def drawRefSignal(self) -> tuple[np.ndarray, np.ndarray]:
        '''
        Generates 1 read signal starting from 3' (RNA) end of the reference and stopping after a uniformly drawn number of bases
        '''
        return self.__drawSignal__()

    def __drawSignal__(self, stop : int = None) -> tuple[np.ndarray, np.ndarray]:
        '''
        Generates a read signal from the whole reference

        Parameters
        ----------
        stop : int
            simulate signal from the 3' end to this base position (base included and 1-based)

        Returns
        -------
        sim_signals : np.ndarray
            a numpy array representing the simulated signal according the given kmer_model
        borders : np.ndarray
            an array containing the segment borders starting with 0
        '''
        
        if stop is not None:
            reference = self.reference[:stop]
        else:
            reference = self.reference

        assert len(reference) > 4, f'Reference sequence too small ({len(reference)}) for simulation'

        # current length of the simulated signal
        signal_pointer = 0
        init_len = len(reference) * int(self.exp_m)
        sim_signal = np.zeros(init_len, dtype = float)

        border_pionter = 1
        borders = np.zeros(len(reference) - 3, dtype = int)

        for n in range(len(reference) - 5):
            kmer = reference[n:n+5]
            segment_length = self.__drawSegmentLength__()

            if signal_pointer + segment_length >= len(sim_signal):
                sim_signal = np.append(sim_signal, np.zeros(init_len, dtype = float))

            sim_signal[signal_pointer : signal_pointer + segment_length] = np.random.normal(self.__drawSegmentSignal__(kmer, segment_length))
            signal_pointer += segment_length

            borders[border_pionter] = signal_pointer
            border_pionter += 1

        return sim_signal[:signal_pointer], borders