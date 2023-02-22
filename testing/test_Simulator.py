#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import sys
sys.path.append("..")
from ..ont_simsig.Simulator import RNASimulator

def testLength():
    rna = RNASimulator(refLength = 1943)
    assert rna.getRefLength() == 1943

def testAlphabet():
    alphabet='ACG'
    rna = RNASimulator(alphabet=alphabet)
    assert set(rna.getReference()) <= set(alphabet)

def testLengthPrefix():
    prefix = 'ACGCG' * 20
    rna = RNASimulator(prefix=prefix)
    assert rna.getRefLength() == 2100

def testLengthSuffix():
    suffix = 'ACGCG' * 20
    rna = RNASimulator(suffix=suffix)
    assert rna.getRefLength() == 2100

def testLengthPrefixSuffix():
    prefix = suffix = 'ACGCG' * 20
    rna = RNASimulator(prefix=prefix, suffix=suffix)
    assert rna.getRefLength() == 2200

def testReferencePrefixSuffix():
    prefix = suffix = 'ACGCG' * 20
    ref = 'CAGCTAGTCGACTAGCTA'
    rna = RNASimulator(reference = ref, prefix=prefix, suffix=suffix)
    assert rna.getReference() == prefix + ref + suffix

def testReference():
    ref = 'CAGCTAGTCGACTAGCTA'
    rna = RNASimulator(reference = ref)
    assert rna.getReference() == ref