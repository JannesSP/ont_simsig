import pandas as pd
import os
from Reference import RNASimulator
from Writer import RNAWriter

# Loading model
kmer_model_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'template_median69pA.model')
df = pd.read_csv(kmer_model_file, sep='\t')
kmer_dict = {key : (mean, std) for key, mean, std in zip(df['kmer'], df['level_mean'], df['level_stdv'])}

print('Simulate RNA reads')
rna = RNASimulator(kmer_dict)
reference = rna.getReference()
signals = rna.drawRefSignals(100)

writer = RNAWriter(reference, path=os.path.join('data','simulation'), dedicated_filename='test')
writer.writeReads(signals)