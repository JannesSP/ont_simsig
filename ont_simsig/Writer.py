# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from datetime import datetime
from os.path import join, exists, splitext, basename
from os import makedirs
from typing import Iterable, Tuple

import h5py
import numpy as np


class RNAWriter():
    '''
    Class to write read signals into the multi FAST5 format
    '''
    
    def __init__(self, reference : str, path : str = '.', dedicatedFilename : str = None, barcoded : bool = False, batchsize : int = 4000, tag : str = '', header : str = None):
        '''
        Parameters
        ----------
        reference : str
            reference string to write into FAST5 file
        dedicated_filename : str
            name of the FAST5 file to write to, without extension!
        barcoded : bool
            changes a flag in to FAST5 files
        batchsize : int
            number of signals to write into the FAST5 file at once
        '''
        self.batch = 0
        self.readNum = 0
        self.startTime = 0
        self.header = header
        self.reference = reference
        self.batchsize = batchsize
        self.barcoded = barcoded
        self.date = datetime.now().strftime("%Y%m%d")
        self.datetimeClean = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.datetime = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")

        path = join(path, f'RNA_simulation_{self.datetimeClean}', tag)
        if dedicatedFilename is not None:
            self.filename = join(path, dedicatedFilename)
        else:
            self.filename = join(path, f'RNA_simulation_{self.datetimeClean}_batch')

        if not exists(path):
            makedirs(path)
        
        self.__initH5()
        self.__initSum()
        self.__writeRefFasta()
        self.__initReadFasta()

    def __initReadFasta(self) -> None:
        self.fasta = open(f'{self.filename}{self.batch}.fasta', 'w')

    def __closeReadFasta(self) -> None:
        self.fasta.close()

    def __initH5(self) -> None:
        self.h5 = h5py.File(f'{self.filename}{self.batch}.fast5', 'w')
        self.h5.attrs.create('file_version', data=np.bytes_('2.2'))
        self.h5.attrs.create('file_type', data=np.bytes_('multi-read'))
        # EXTRA INFORMATION
        self.h5.attrs.create('num_nucleotides', data=len(self.reference), dtype=np.uint16)
        self.h5.attrs.create('Reference', data=np.string_(self.reference))
        
    def __initSum(self) -> None:
        self.sum = open(f'{self.filename}_sequencing_summary.txt', 'w')
        self.sum.write('filename_fastq\tfilename_fast5\tread_id\trun_id\tchannel\tmux\tstart_time\tduration\tpore_type\texperiment_id\tsample_id\tend_reason\n')

    def __writeRefFasta(self) -> None:
        with open(f'{splitext(self.filename)[0]}_reference.fasta', 'w') as f:
            if self.header is not None:
                f.write(f'>{self.header}\n{self.reference}\n')
            else:
                f.write(f'>{basename(self.filename)}\n{self.reference}\n')

    def getFilename(self) -> str:
        return self.filename

    def writeRead(self, simSignal : Tuple[np.ndarray, np.ndarray]) -> None:
        '''
        Write one signal to the FAST5 file
        Parameters
        ----------
        signal : tuple[np.ndarray, np.ndarray]
            pA signal to write into the FAST5 file
            borders of the segments to write into the FAST5 file
        '''
        self.writeReads([simSignal])
    
    def writeReads(self, simSignals : Iterable[Tuple[np.ndarray, np.ndarray, list]]) -> None:
        '''
        Write multiple signals to the FAST5 file
        Parameters
        ----------
        simSignal : Iterable[tuple[np.ndarray, np.ndarray]]
            pA signals to write into the FAST5 file
            borders of the segments to write into the FAST5 file
        '''
        for num, (signal, borders, read) in enumerate(simSignals):
            
            if (num+1)%100==0:
                print(f'Writing read {self.readNum + 1}\{len(simSignals)} in batch {self.batch} ...', end = '\r')
            
            if self.readNum%self.batchsize == 0 and self.readNum != 0:
                self.h5.close()
                self.batch += 1
                self.h5 = h5py.File(f'{self.filename}{self.batch}.fast5', 'w')
                self.h5.attrs.create("file_version", data=np.bytes_('2.2'))
                self.h5.attrs.create("file_type", data=np.bytes_('multi-read'))
                # EXTRA INFORMATION
                self.h5.attrs.create('num_nucleotides', data=len(self.reference), dtype=np.uint16)
                self.h5.attrs.create("Reference", data=np.string_(self.reference))

                self.__closeReadFasta()
                self.__initReadFasta()
            
            readId = 'sim-' + str(self.readNum)
            fast5Id = 'read_' + readId
            channelNumber = str(np.random.randint(1, 513))

            self.fasta.write(readId + '\n' + ''.join(read[::-1]) + '\n')

            read = self.h5.create_group(fast5Id)
            read.attrs.create('pore_type', data=np.bytes_('not_set'))
            read.attrs.create('run_id', data=np.bytes_('rna_simulation'))
            
            raw = read.create_group('Raw')
            raw.attrs.create('duration', data=len(signal), dtype=np.uint32)
            raw.attrs.create('end_reason', data=5, dtype=np.uint8)
            raw.attrs.create('median_before', data=np.random.normal(217.59, 22.53), dtype=np.float64) # approximated from some real data
            raw.attrs.create('read_id', data=np.bytes_(readId))
            raw.attrs.create('read_number', data=self.readNum, dtype=np.int32)
            raw.attrs.create('start_mux', data=0, dtype=np.uint8)
            raw.attrs.create('start_time', data=self.startTime, dtype=np.uint64)
            raw.create_dataset('Signal', data=np.ceil(signal * (8192/1119.071533203125) + 0), dtype=np.int16) # from sarscov2 kiel data 22195
            # ATTENTION: only looked into 1 fast5 batch file, 8192 is always the same, the range of 1119.... could change, I just took one value from it
            # offset is set to 0, it changes by small value of [-10, 10] as far as I saw in the data
            
            # EXTRA INFORMATION
            raw.attrs.create('num_segments', data=len(borders) - 1, dtype=np.uint64)
            raw.create_dataset('Borders', data=borders, dtype=np.uint64)
            
            channelId = read.create_group('channel_id')
            channelId.attrs.create('channel_number', data=np.bytes_(channelNumber))
            # current = (Dacs + offset ) * range / digitisation = (Dacs + 0) * (1 / 1) <=> current = Dacs in this case
            channelId.attrs.create('digitisation', data=8192, dtype=np.float64) # from sarscov2 kiel data 22195
            channelId.attrs.create('offset', data=0, dtype=np.float64) # from sarscov2 kiel data 22195 (can change, not always 0)
            channelId.attrs.create('range', data=1119.071533203125, dtype=np.float64) # from sarscov2 kiel data 22195
            channelId.attrs.create('sampling_rate', data=3012, dtype=np.float64) # always set to 3012 Hz
            
            contextTags = read.create_group('context_tags')
            contextTags.attrs.create('barcoding_enabled', data=np.bytes_(self.barcoded))
            contextTags.attrs.create('experiment_duration_set', data=np.bytes_('4320'))
            contextTags.attrs.create('experiment_type', data=np.bytes_('rna')) # writer only used for RNA
            contextTags.attrs.create('local_basecalling', data=np.bytes_(False))
            contextTags.attrs.create('package', data=np.bytes_('bream4'))
            contextTags.attrs.create('package_version', data=np.bytes_('6.1.4'))
            contextTags.attrs.create('sample_frequency', data=np.bytes_('3012'))
            contextTags.attrs.create('sequencing_kit', data=np.bytes_('sqk-rna002'))
            
            # DATA from sarscov_kiel 22195 sample
            trackingId = read.create_group('tracking_id')
            trackingId.attrs.create('asic_id', data=np.bytes_('614860902'))
            trackingId.attrs.create('asic_id_eeprom', data=np.bytes_('5532807'))
            trackingId.attrs.create('asic_temp', data=np.bytes_('23.729523'))
            trackingId.attrs.create('asic_version', data=np.bytes_('IA02D'))
            trackingId.attrs.create('auto_update', data=np.bytes_(False))
            trackingId.attrs.create('auto_update_source', data=np.bytes_('https://mirror.oxfordnanoportal.com/software/MinKNOW/'))
            trackingId.attrs.create('bream_is_standard', data=np.bytes_(False))
            trackingId.attrs.create('configuration_version', data=np.bytes_('4.1.15'))
            trackingId.attrs.create('device_id', data=np.bytes_('MN20569'))
            trackingId.attrs.create('device_type', data=np.bytes_('minion'))
            trackingId.attrs.create('distribution_status', data=np.bytes_('stable'))
            trackingId.attrs.create('distribution_version', data=np.bytes_('20.10.3'))
            trackingId.attrs.create('exp_script_name', data=np.bytes_('sequencing/sequencing_MIN106_RNA:FLO-MIN106:SQK-RNA002'))
            trackingId.attrs.create('exp_script_purpose', data=np.bytes_('sequencing_run'))
            trackingId.attrs.create('exp_start_time', data=np.bytes_(self.datetime))
            trackingId.attrs.create('flow_cell_id', data=np.bytes_('FAO86549'))
            trackingId.attrs.create('flow_cell_product_code', data=np.bytes_('FLO-MIN106'))
            trackingId.attrs.create('guppy_version', data=np.bytes_('4.2.2+effbaf8'))
            trackingId.attrs.create('heatsink_temp', data=np.bytes_('34.000000'))
            trackingId.attrs.create('hostname', data=np.bytes_('simulation_pc'))
            trackingId.attrs.create('installation_type', data=np.bytes_('nc'))
            trackingId.attrs.create('local_firmware_file', data=np.bytes_('1'))
            trackingId.attrs.create('operating_system', data=np.bytes_('simulation_os'))
            trackingId.attrs.create('protocol_group_id', data=np.bytes_(f'{self.date}_simulation_group_id'))
            trackingId.attrs.create('protocol_run_id', data=np.bytes_('62ce610f-8b14-4236-a120-7d940187a02d'))
            trackingId.attrs.create('protocols_version', data=np.bytes_('6.1.4'))
            trackingId.attrs.create('run_id', data=np.bytes_(f'{self.date}_simulation_run_id'))
            trackingId.attrs.create('sample_id', data=np.bytes_(f'{self.date}_simulation_sample_id'))
            trackingId.attrs.create('usb_config', data=np.bytes_('MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto'))
            trackingId.attrs.create('version', data=np.bytes_('4.1.2'))
            
            self.sum.write(f'-\t{self.filename}{self.batch}.fast5\t{readId}\t{self.date}_simulation_run_id\t{channelNumber}\t0\t{self.startTime}\t{len(signal)}\tnot_set\t{self.date}_simulation_run_id\t{self.date}_simulation_sample_id\tsignal_positive\n')
            
            self.startTime += len(signal)
            self.readNum += 1
        
        self.__closeReadFasta()
        self.h5.close()
        print(f'\nDone writing {len(simSignals)} reads')
            
    def close(self):
        '''
        Close FAST5 file
        Close summary file
        '''
        self.h5.close()
        self.sum.close()
