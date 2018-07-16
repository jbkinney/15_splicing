#!/usr/bin/env python

import pandas as pd
import numpy as np

import os
import sys
import re
import pdb
import time
from pyfasta import Fasta

sys.path.append('src')
from utils import make_iupac_seqs

print('Running analyze_kmers.py. YOU SHOULD NOT HAVE TO RUN THIS.'
      ' RESULTS ARE ALREADY IN DATA DIRECTORY.')

print('Tabulating kmers in genome of the form NNNGYNNNN')

# Set input and output files
hg38_fasta_file = 'data/hg38.fa'
out_file = 'data/kmer_counts.txt'

# Load human fasta
hg38 = Fasta(hg38_fasta_file)

# Function to compute the reverse complement of a sequence
complement_dict = ''.maketrans('ATCGN','TAGCN')
def rc(seq):
    s = seq.upper().translate(complement_dict)
    return s[::-1]

# Get list of chromosomes
chroms = [l for l in hg38.keys() if ('_' not in l) and ('M' not in l)]
chroms.sort()

# Build empty dict for storing kmers
kmers = make_iupac_seqs('NNNGYNNNN')
kmer_counts_dict = dict((kmer,0) for kmer in kmers)

# For each chromosome, tally k-mers
t = time.time()
for chrom in chroms:
    print('Processing %s...'%chrom)
    seq = hg38[chrom][:]
    
    # Tally kmers in forward direction
    print('forward...')
    fwd_iterator = re.finditer('[ACGTacgt]{3}G[CT][ACGTacgt]{4}', seq[:])
    for match in fwd_iterator:
        kmer = match.group().upper()
        kmer_counts_dict[kmer] += 1
        
    # Tally kmers in reverse direction
    print('reverse...')
    rev_iterator = re.finditer('[ACGTacgt]{4}[AG]C[ACGTacgt]{3}', seq[:])
    for match in rev_iterator:
        kmer = rc(match.group().upper())
        kmer_counts_dict[kmer] += 1
    
kmers = list(kmer_counts_dict.keys())
counts = list(kmer_counts_dict.values())

df = pd.DataFrame()
df['kmer'] = kmers
df['count'] = counts

df.to_csv(out_file, sep='\t', index=False)
print('Done!')
