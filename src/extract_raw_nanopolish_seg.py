# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import sys
import copy
import h5py
import os
import numpy as np
# TODO Logger is another project accessible with export PYTHONPATH='~/projects'
from Logger.Logger import Logger

logger = Logger()

def str2int(c):
    alph = ['A', 'C', 'G', 'T']
    try:
        return alph.index(c)
    except ValueError:
        return -1

def read_line(line, header = False):
    line_dic = {}

    if header:

        line_dic['contig'] = ''
        line_dic['position'] = -1
        line_dic['ref_kmer'] = ''
        line_dic['read_index'] = -1
        line_dic['strand'] = ''
        line_dic['event_index'] = -1
        line_dic['event_mean'] = 0.0
        line_dic['event_std'] = 0.0
        line_dic['event_len'] = 0.0
        line_dic['model_kmer'] = ''
        line_dic['model_mean'] = 0.0
        line_dic['model_std'] = 0.0
        line_dic['standardized_level'] = 0
        line_dic['start_idx'] = -1
        line_dic['end_idx'] = -1

    else:

        line = line.strip().split()

        line_dic['contig'] = line[0]
        line_dic['position'] = int(line[1])
        line_dic['ref_kmer'] = line[2]
        line_dic['read_index'] = int(line[3])
        line_dic['strand'] = line[4]
        line_dic['event_index'] = int(line[5])
        line_dic['event_mean'] = float(line[6])
        line_dic['event_std'] = float(line[7])
        line_dic['event_len'] = float(line[8])
        line_dic['model_kmer'] = line[9]
        line_dic['model_mean'] = float(line[10])
        line_dic['model_std'] = float(line[11])
        line_dic['standardized_level'] = float(line[12])
        line_dic['start_idx'] = int(line[13])
        line_dic['end_idx'] = int(line[14])

    return line_dic

# PARAMS
parser = ArgumentParser(
    formatter_class=ArgumentDefaultsHelpFormatter,
    add_help=False
)

parser.add_argument("nanopolish_summary")
parser.add_argument("nanopolish_result")
parser.add_argument("output_hdf5", help='output filename')

args = parser.parse_args()

nano_sum = args.nanopolish_summary
nano_res = args.nanopolish_result
output_hdf5 = args.output_hdf5

assert nano_sum.endswith('.csv') and os.path.exists(nano_sum)
assert nano_res.endswith('.csv') and os.path.exists(nano_res)
assert output_hdf5.endswith('.hdf5'), 'Please set a hdf5 file as output.'
assert not os.path.exists(output_hdf5), f'{output_hdf5} already exists!'

#### read summary file with read_indices and read_ids
read_index2ID = {}
logger.printLog(f'Start loading ids from {nano_sum} ...')
with open(nano_sum, 'r') as sum:

    # skip header
    next(sum)

    for line in sum:

        r_index, r_ID = line.strip().split()[:2]

        read_index2ID[int(r_index)] = r_ID
logger.printLog('Done')

output_fh = h5py.File(output_hdf5, 'w', libver='latest')

logger.printLog(f'Start preparing data with nanopolish segmentation from {nano_res} ...')
#### read result with and get segmentation
with open(nano_res, 'r') as res:

    # read first line
    cur_event = read_line(res.readline(), True)
    segmentation = []

    for l_idx, line in enumerate(res):

        if (l_idx + 1) % 100000 == 0:
            sys.stderr.write(f'\rLine {l_idx + 1}')

        line = read_line(line)

        # new read
        if cur_event['read_index'] != line['read_index'] and cur_event['read_index'] != -1:

            seg_array = np.array(segmentation)

            seg_array = seg_array[seg_array[:, 0].argsort()]
            
            if read_index2ID[cur_event['read_index']] not in output_fh:
            
                output_fh[read_index2ID[cur_event['read_index']]] = seg_array[np.unique(seg_array[:, 0], axis = 0, return_index=True)[1]]

            else:
            
                temp_seg = output_fh[read_index2ID[cur_event['read_index']]]
                seg_array = np.append(seg_array, temp_seg, axis = 0)
                seg_array = seg_array[seg_array[:, 0].argsort()]
                del output_fh[read_index2ID[cur_event['read_index']]]

                unique_seg_array = []
                prev_seg = seg_array[0, 0]
                base_list = [seg_array[0, 1]]

                # if the segmentation is not unique, we need to merge them, take most abundant base
                for seg, base in seg_array[1:]:

                    if prev_seg != seg:

                        if max(base_list, key = base_list.count) == -1 and len(list(filter((-1).__ne__, base_list))) > 0:

                            base_list = list(filter((-1).__ne__, base_list))

                        unique_seg_array.append([prev_seg, max(base_list, key = base_list.count)])
                        
                        base_list = []

                    base_list.append(base)
                    prev_seg = seg

                if max(base_list, key = base_list.count) == -1 and len(list(filter((-1).__ne__, base_list))) > 0:

                    base_list = list(filter((-1).__ne__, base_list))

                unique_seg_array.append([prev_seg, max(base_list, key = base_list.count)])
                unique_seg_array = np.array(unique_seg_array)

                if not len(unique_seg_array) == len(np.unique(unique_seg_array[:, 0], axis = 0)):
                    print(unique_seg_array)
                    assert len(unique_seg_array) == len(np.unique(unique_seg_array[:, 0], axis = 0))

                output_fh[read_index2ID[cur_event['read_index']]] = unique_seg_array[np.unique(unique_seg_array[:, 0], axis = 0, return_index=True)[1]]

            segmentation = []
            cur_event = copy.deepcopy(line)

        # new position in reference
        elif cur_event['position'] != line['position'] and cur_event['read_index'] != -1:

            segmentation.append([cur_event['end_idx'], str2int(cur_event['model_kmer'][2])])

            cur_event = copy.deepcopy(line)

        else:

            # concatenate events
            temp_start = cur_event['start_idx']
            temp_end = cur_event['end_idx']
            cur_event = copy.deepcopy(line)

            if temp_start == line['end_idx']:

                cur_event['end_idx'] = temp_end

            elif temp_end == line['start_idx']:

                cur_event['start_idx'] = temp_start