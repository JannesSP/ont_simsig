from typing import Iterable
import h5py
import random
from datetime import datetime
import numpy as np

class RNAWriter():
    '''
    Class to write read signals into the multi FAST5 format
    '''
    
    def __init__(self, barcoded : bool = False, batchsize : int = 4000):
        '''
        Parameters
        ----------
        barcoded : bool
            changes a flag in to FAST5 files
        batchsize : int
            number of signals to write into the FAST5 file at once
        '''
        
        self.batch = 0
        self.read_num = 0
        self.start_time = 0
        
        self.batchsize = batchsize
        self.barcoded = barcoded
        
        self.date = datetime.now().strftime("%Y%m%d")
        self.datetime = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
        self.h5 = h5py.File(f'RNA_simulation_{self.date}_batch{self.batch}.fast5', 'w')

    def writeReas(self, signal : np.ndarray) -> None:
        '''
        Write one signal to the FAST5 file
        Parameters
        ----------
        signal : np.ndarray
            pA signal to write into the FAST5 file
        '''
        self.writeReads([signal])    
    
    def writeReads(self, signals : Iterable[np.ndarray]) -> None:
        '''
        Write multiple signals to the FAST5 file
        Parameters
        ----------
        signals : Iterable[np.ndarray]
            pA signals to write into the FAST5 file
        '''
        for signal in signals:
            
            if (self.read_num + 1)%self.batchsize == 0:
                self.h5.close()
                self.batch += 1
                self.h5 = h5py.File(f'RNA_simulation_{self.date}_batch{self.batch}.fast5', 'w')
            
            readid = 'read_' + str(self.read_num)
            channel_number = random.randint(1, 512)
            
            self.h5.attrs.create("file_version", data=b'2.2', dtype=np.bytes_)
            self.h5.attrs.create("file_type", data=b'multi-read', dtype=np.bytes_)

            read = self.h5.create_group(readid)
            read.attrs.create('pore_type', data=b'not_set', dtype=np.bytes_)
            read.attrs.create('run_id', data=b'rna_simulation', dtype=np.bytes_)
            
            raw = read.create_group('Raw')
            raw.attrs.create('duration', data=len(signal), dtype=np.uint32)
            raw.attrs.create('end_reason', data=5, dtype=np.uint8)
            raw.attrs.create('median_before', data=np.random.normal(217.59, 22.53), dtype=np.float64) # approximated from some real data
            raw.attrs.create('read_id', data=self.read_num, dtype=np.bytes_)
            raw.attrs.create('read_number', data=self.read_num, dtype=np.int32)
            raw.attrs.create('start_mux', data=0, dtype=np.uint8)
            raw.attrs.create('start_time', data=self.start_time, dtype=np.uint64)
            raw.create_dataset('Signal', data=signal, dtype=np.float64)
            
            channel_id = read.create_group('channel_id')
            channel_id.attrs.create('channel_number', data=channel_number, dtype=np.bytes_)
            # current = (Dacs + offset ) * range / digitisation = (Dacs + 0) * (1 / 1) <=> current = Dacs in this case
            channel_id.attrs.create('digitisation', data=1, dtype=np.float64) # always set to 1
            channel_id.attrs.create('offset', data=0, dtype=np.float64) # always set to 0
            channel_id.attrs.create('range', data=1, dtype=np.float64) # always set to 1
            channel_id.attrs.create('sampling_rate', data=3012, dtype=np.float64) # always set to 3012 Hz
            
            context_tags = read.create_group('context_tags')
            context_tags.attrs.create('barcoding_enabled', data=self.barcoded, dtype=np.bytes_)
            context_tags.attrs.create('experiment_duration_set', data=b'4320', dtype=np.bytes_)
            context_tags.attrs.create('experiment_type', data=b'rna', dtype=np.bytes_) # writer only used for RNA
            context_tags.attrs.create('local_basecalling', data=False, dtype=np.bytes_)
            context_tags.attrs.create('package', data=b'bream4', dtype=np.bytes_)
            context_tags.attrs.create('package_version', data=b'6.1.4', dtype=np.bytes_)
            context_tags.attrs.create('sample_frequency', data=b'3012', dtype=np.bytes_)
            context_tags.attrs.create('sequencing_kit', data=b'sqk-rna002', dtype=np.bytes_)
            
            tracking_id = read.create_group('tracking_id')
            tracking_id.attrs.create('asic_id', data=b'614860902', dtype=np.bytes_)
            tracking_id.attrs.create('asic_id_eeprom', data=b'5532807', dtype=np.bytes_)
            tracking_id.attrs.create('asic_temp', data='23.729523', dtype=np.bytes_)
            tracking_id.attrs.create('asic_version', data=b'IA02D', dtype=np.bytes_)
            tracking_id.attrs.create('auto_update', data=False, dtype=np.bytes_)
            tracking_id.attrs.create('auto_update_source', data=b'https://mirror.oxfordnanoportal.com/software/MinKNOW/', dtype=np.bytes_)
            tracking_id.attrs.create('bream_is_standard', data=False, dtype=np.bytes_)
            tracking_id.attrs.create('configuration_version', data=b'4.1.15', dtype=np.bytes_)
            tracking_id.attrs.create('device_id', data=b'MN20569', dtype=np.bytes_)
            tracking_id.attrs.create('device_type', data=b'minion', dtype=np.bytes_)
            tracking_id.attrs.create('distribution_status', data=b'stable', dtype=np.bytes_)
            tracking_id.attrs.create('distribution_version', data=b'20.10.3', dtype=np.bytes_)
            tracking_id.attrs.create('exp_script_name', data=b'sequencing/sequencing_MIN106_RNA:FLO-MIN106:SQK-RNA002', dtype=np.bytes_)
            tracking_id.attrs.create('exp_script_purpose', data=b'sequencing_run', dtype=np.bytes_)
            tracking_id.attrs.create('exp_start_time', data=self.datetime, dtype=np.bytes_)
            tracking_id.attrs.create('flow_cell_id', data=b'FAO86549', dtype=np.bytes_)
            tracking_id.attrs.create('flow_cell_product_code', data=b'FLO-MIN106', dtype=np.bytes_)
            tracking_id.attrs.create('guppy_version', data=b'4.2.2+effbaf8', dtype=np.bytes_)
            tracking_id.attrs.create('heatsink_temp', data=b'34.000000', dtype=np.bytes_)
            tracking_id.attrs.create('hostname', data=b'simulation_pc', dtype=np.bytes_)
            tracking_id.attrs.create('installation_type', data=b'nc', dtype=np.bytes_)
            tracking_id.attrs.create('local_firmware_file', data=b'1', dtype=np.bytes_)
            tracking_id.attrs.create('operating_system', data=b'simulation_os', dtype=np.bytes_)
            tracking_id.attrs.create('protocol_group_id', data=f'{self.date}_simulation_group_id', dtype=np.bytes_)
            tracking_id.attrs.create('protocol_run_id', data=b'62ce610f-8b14-4236-a120-7d940187a02d', dtype=np.bytes_)
            tracking_id.attrs.create('protocols_version', data=b'6.1.4', dtype=np.bytes_)
            tracking_id.attrs.create('run_id', data=f'{self.date}_simulation_run_id', dtype=np.bytes_)
            tracking_id.attrs.create('sample_id', data=f'{self.date}_simulation_sample_id', dtype=np.bytes_)
            tracking_id.attrs.create('usb_config', data=b'MinION_fx3_1.1.1_ONT#MinION_fpga_1.1.0#bulk#Auto', dtype=np.bytes_)
            tracking_id.attrs.create('version', data=b'4.1.2', dtype=np.bytes_)
            
            self.start_time += len(signal)
            self.read_num += 1
            
    def closeWrite(self):
        '''
        Close FAST5 file
        '''
        self.h5.close()