from src.Reference import RNASimulator
import os
import pandas as pd
import re

def readData():
    kmer_model_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'template_median69pA.model')
    df = pd.read_csv(kmer_model_file, sep='\t')
    return {key : (mean, std) for key, mean, std in zip(df['kmer'], df['level_mean'], df['level_stdv'])}

kmer_dict = readData()

def testLength():
    rna = RNASimulator(kmer_dict, length = 1943)
    assert rna.getRefLength() == 1943

def testAlphabet():
    alphabet='ACG'
    rna = RNASimulator(kmer_dict, alphabet=alphabet)
    assert set(rna.getReference()) <= set(alphabet)

def testLengthPrefix():
    prefix = 'ACGCG' * 20
    rna = RNASimulator(kmer_dict, prefix=prefix)
    assert rna.getRefLength() == 2100

def testLengthSuffix():
    suffix = 'ACGCG' * 20
    rna = RNASimulator(kmer_dict, suffix=suffix)
    assert rna.getRefLength() == 2100

def testLengthPrefixSuffix():
    prefix = suffix = 'ACGCG' * 20
    rna = RNASimulator(kmer_dict, prefix=prefix, suffix=suffix)
    assert rna.getRefLength() == 2200

def testReferencePrefixSuffix():
    prefix = suffix = 'ACGCG' * 20
    ref = 'CAGCTAGTCGACTAGCTA'
    rna = RNASimulator(kmer_dict, reference = ref, prefix=prefix, suffix=suffix)
    assert rna.getReference() == prefix + ref + suffix

def testReference():
    ref = 'CAGCTAGTCGACTAGCTA'
    rna = RNASimulator(kmer_dict, reference = ref)
    assert rna.getReference() == ref