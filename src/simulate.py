import pandas as pd
# import numpy as np
import matplotlib.pyplot as plt
import os
# import random
# from scipy import stats
from Reference import RNASimulator

kmer_model_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'template_median69pA.model')
df = pd.read_csv(kmer_model_file, sep='\t')
kmer_dict = {key : (mean, std) for key, mean, std in zip(df['kmer'], df['level_mean'], df['level_stdv'])}

rna = RNASimulator(kmer_dict)

test_signal, test_borders = rna.drawRefSignal()

plt.plot(test_signal)
plt.show()