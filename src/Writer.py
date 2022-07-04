from datetime import datetime
from os.path import join, exists
from os import makedirs
from typing import Iterable

import h5py
import numpy as np


class RNAWriter():
    '''
    Class to write read signals into the multi FAST5 format
    '''
    
    def __init__(self, reference : str, path : str = '.', dedicated_filename : str = None, barcoded : bool = False, batchsize : int = 4000):
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
        self.read_num = 0
        self.start_time = 0
        self.reference = reference
        self.batchsize = batchsize
        self.barcoded = barcoded
        self.date = datetime.now().strftime("%Y%m%d")
        self.datetime_clean = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.datetime = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")

        path = join(path, f'RNA_simulation_{self.datetime_clean}')
        if dedicated_filename is not None:
            self.filename = join(path, dedicated_filename)
        else:
            self.filename = join(path, f'RNA_simulation_{self.datetime_clean}_batch')

        if not exists(path):
            makedirs(path)
        
        self.__initH5()
        self.__initSum()
        self.__writeRefFasta()

    def __initH5(self) -> None:
        self.h5 = h5py.File(f'{self.filename}{self.batch}.fast5', 'w')
        self.h5.attrs.create('file_version', data=np.bytes_('2.2'))
        self.h5.attrs.create('file_type', data=np.bytes_('multi-read'))
        # EXTRA INFORMATION
        self.h5.attrs.create('num_nucleotides', data=len(self.reference), dtype=np.uint16)
        self.h5.create_dataset('Reference', data=np.string_(self.reference))
        
    def __initSum(self) -> None:
        self.sum = open(f'{self.filename}_sequencing_summary.txt', 'w')
        self.sum.write('filename_fastq\tfilename_fast5\tread_id\trun_id\tchannel\tmux\tstart_time\tduration\tpore_type\texperiment_id\tsample_id\tend_reason\n')

    def __writeRefFasta(self) -> None:
        with open(f'{self.filename}_reference.fasta', 'w') as f:
            f.write(f'>{self.filename}\n{self.reference}\n')

    def getFilename(self) -> str:
        return self.filename

    def writeRead(self, simSignal : tuple[np.ndarray, np.ndarray]) -> None:
        '''
        Write one signal to the FAST5 file
        Parameters
        ----------
        signal : tuple[np.ndarray, np.ndarray]
            pA signal to write into the FAST5 file
            borders of the segments to write into the FAST5 file
        '''
        self.writeReads([simSignal])
    
    def writeReads(self, simSignals : Iterable[tuple[np.ndarray, np.ndarray]]) -> None:
        '''
        Write multiple signals to the FAST5 file
        Parameters
        ----------
        simSignal : Iterable[tuple[np.ndarray, np.ndarray]]
            pA signals to write into the FAST5 file
            borders of the segments to write into the FAST5 file
        '''
        for num, (signal, borders) in enumerate(simSignals):
            
            if (num+1)%10==0:
                print(f'Writing read {self.read_num + 1}\{len(simSignals)} in batch {self.batch} ...', end = '\r')
            
            if self.read_num%self.batchsize == 0 and self.read_num != 0:
                self.h5.close()
                self.batch += 1
                self.h5 = h5py.File(f'{self.filename}{self.batch}.fast5', 'w')
                self.h5.attrs.create("file_version", data=np.bytes_('2.2'))
                self.h5.attrs.create("file_type", data=np.bytes_('multi-read'))
                # EXTRA INFORMATION
                self.h5.attrs.create('num_nucleotides', data=len(self.reference), dtype=np.uint16)
                self.h5.create_dataset("Reference", data=np.string_(self.reference))
            
            readid = 'read_' + str(self.read_num)
            channel_number = str(np.random.randint(1, 513))

            read = self.h5.create_group(readid)
            read.attrs.create('pore_type', data=np.bytes_('not_set'))
            read.attrs.create('run_id', data=np.bytes_('rna_simulation'))
            
            raw = read.create_group('Raw')
            raw.attrs.create('duration', data=len(signal), dtype=np.uint32)
            raw.attrs.create('end_reason', data=5, dtype=np.uint8)
            raw.attrs.create('median_before', data=np.random.normal(217.59, 22.53), dtype=np.float64) # approximated from some real data
            raw.attrs.create('read_id', data=str(self.read_num))
            raw.attrs.create('read_number', data=self.read_num, dtype=np.int32)
            raw.attrs.create('start_mux', data=0, dtype=np.uint8)
            raw.attrs.create('start_time', data=self.start_time, dtype=np.uint64)
            raw.create_dataset('Signal', data=np.ceil(signal * (8192/1119.071533203125) + 0), dtype=np.int16) # from sarscov2 kiel data 22195
            
            # EXTRA INFORMATION
            raw.attrs.create('num_segments', data=len(borders) - 1, dtype=np.uint64)
            raw.create_dataset('Borders', data=borders, dtype=np.uint64)
            
            channel_id = read.create_group('channel_id')
            channel_id.attrs.create('channel_number', data=np.bytes_(channel_number))
            # current = (Dacs + offset ) * range / digitisation = (Dacs + 0) * (1 / 1) <=> current = Dacs in this case
            channel_id.attrs.create('digitisation', data=8192, dtype=np.float64) # from sarscov2 kiel data 22195
            channel_id.attrs.create('offset', data=0, dtype=np.float64) # from sarscov2 kiel data 22195 (can change, not always 0)
            channel_id.attrs.create('range', data=1119.071533203125, dtype=np.float64) # from sarscov2 kiel data 22195
            channel_id.attrs.create('sampling_rate', data=3012, dtype=np.float64) # always set to 3012 Hz
            
            context_tags = read.create_group('context_tags')
            context_tags.attrs.create('barcoding_enabled', data=np.bytes_(self.barcoded))
            context_tags.attrs.create('experiment_duration_set', data=np.bytes_('4320'))
            context_tags.attrs.create('experiment_type', data=np.bytes_('rna')) # writer only used for RNA
            context_tags.attrs.create('local_basecalling', data=np.bytes_(False))
            context_tags.attrs.create('package', data=np.bytes_('bream4'))
            context_tags.attrs.create('package_version', data=np.bytes_('6.1.4'))
            context_tags.attrs.create('sample_frequency', data=np.bytes_('3012'))
            context_tags.attrs.create('sequencing_kit', data=np.bytes_('sqk-rna002'))
            
            tracking_id = read.create_group('tracking_id')
            tracking_id.attrs.create('asic_id', data=np.bytes_('614860902'))
            tracking_id.attrs.create('asic_id_eeprom', data=np.bytes_('5532807'))
            tracking_id.attrs.create('asic_temp', data=np.bytes_('23.729523'))
            tracking_id.attrs.create('asic_version', data=np.bytes_('IA02D'))
            tracking_id.attrs.create('auto_update', data=np.bytes_(False))
            tracking_id.attrs.create('auto_update_source', data=np.bytes_('https://mirror.oxfordnanoportal.com/software/MinKNOW/'))
            tracking_id.attrs.create('bream_is_standard', data=np.bytes_(False))
            tracking_id.attrs.create('configuration_version', data=np.bytes_('4.1.15'))
            tracking_id.attrs.create('device_id', data=np.bytes_('MN20569'))
            tracking_id.attrs.create('device_type', data=np.bytes_('minion'))
            tracking_id.attrs.create('distribution_status', data=np.bytes_('stable'))
            tracking_id.attrs.create('distribution_version', data=np.bytes_('20.10.3'))
            tracking_id.attrs.create('exp_script_name', data=np.bytes_('sequencing/sequencing_MIN106_RNA:FLO-MIN106:SQK-RNA002'))
            tracking_id.attrs.create('exp_script_purpose', data=np.bytes_('sequencing_run'))
            tracking_id.attrs.create('exp_start_time', data=np.bytes_(self.datetime))
            tracking_id.attrs.create('flow_cell_id', data=np.bytes_('FAO86549'))
            tracking_id.attrs.create('flow_cell_product_code', data=np.bytes_('FLO-MIN106'))
            tracking_id.attrs.create('guppy_version', data=np.bytes_('4.2.2+effbaf8'))
            tracking_id.attrs.create('heatsink_temp', data=np.bytes_('34.000000'))
            tracking_id.attrs.create('hostname', data=np.bytes_('simulation_pc'))
            tracking_id.attrs.create('installation_type', data=np.bytes_('nc'))
            tracking_id.attrs.create('local_firmware_file', data=np.bytes_('1'))
            tracking_id.attrs.create('operating_system', data=np.bytes_('simulation_os'))
            tracking_id.attrs.create('protocol_group_id', data=np.bytes_(f'{self.date}_simulation_group_id'))
            tracking_id.attrs.create('protocol_run_id', data=np.bytes_('62ce610f-8b14-4236-a120-7d940187a02d'))
            tracking_id.attrs.create('protocols_version', data=np.bytes_('6.1.4'))
            tracking_id.attrs.create('run_id', data=np.bytes_(f'{self.date}_simulation_run_id'))
            tracking_id.attrs.create('sample_id', data=np.bytes_(f'{self.date}_simulation_sample_id'))
            tracking_id.attrs.create('usb_config', data=np.bytes_('MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto'))
            tracking_id.attrs.create('version', data=np.bytes_('4.1.2'))
            
            self.sum.write(f'-\t{self.filename}{self.batch}.fast5\t{self.read_num}\t{self.date}_simulation_run_id\t{channel_number}\t0\t{self.start_time}\t{len(signal)}\tnot_set\t{self.date}_simulation_run_id\t{self.date}_simulation_sample_id\tsignal_positive\n')
            
            self.start_time += len(signal)
            self.read_num += 1
        
        print(f'\nDone writing {len(simSignals)} reads')
            
    def close(self):
        '''
        Close FAST5 file
        Close summary file
        '''
        self.h5.close()
        self.sum.close()
