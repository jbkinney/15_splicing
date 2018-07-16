#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import os
import sys
import re
import pdb
import glob
import time
from pyfasta import Fasta

print('Running analyze_exac_variants.py...')

# Function to compute the reverse complement of a sequence
complement = ''.maketrans('ATCGN', 'TAGCN')
def rc(seq):
    return seq.upper().translate(complement)[::-1]

# Specify input file names
hg_fasta_file = 'data/hg19.fa'
hg_exon_file = 'data/hg19_exons.txt'
snp_file = 'data/exac.splice.filtered.txt'

# Specify output files
exac_snp_file = 'output/exac_5ss_snps.txt'

# Specify chromosomes to which we'll restrict analysis
chromosomes = ['chr%d'%n for n in range(1,23)]+['chrX','chrY']

### Load exon coordinates
col_names = ['chromosome','start','stop','id','strand']
col_indices = [0,1,2,3,5]
exon_coords_df = pd.read_csv(hg_exon_file,sep='\t',names=col_names,usecols=col_indices)
exon_coords_df = exon_coords_df[exon_coords_df['chromosome'].isin(chromosomes)]
exon_coords_df = exon_coords_df[exon_coords_df['strand'].isin(['+','-'])]
print('Loaded position data for %d exons.'%len(exon_coords_df))

# Compute the coordinates of the (potential) 5' ss at the end of each exon
ss_df = exon_coords_df[['chromosome','strand']].copy()
ss_df['start'] = 0
ss_df['stop'] = 0

# Compute ss coordinates on plus strand
plus_indices = exon_coords_df['strand']=='+'
ss_df.loc[plus_indices,'start'] = exon_coords_df.loc[plus_indices,'stop']-3
ss_df.loc[plus_indices,'stop'] = exon_coords_df.loc[plus_indices,'stop']+6

# Compute ss coordinates on minus strand
minus_indices = exon_coords_df['strand']=='-'
ss_df.loc[minus_indices,'start'] = exon_coords_df.loc[minus_indices,'start']-6
ss_df.loc[minus_indices,'stop'] = exon_coords_df.loc[minus_indices,'start']+3

# Drop duplicates splice sites
ss_df.drop_duplicates(inplace=True)
ss_df.reset_index(drop=True,inplace=True)

# Load SNP data
snp_df = pd.read_csv(snp_file, 
                     usecols=['chrom','pos','ref','alt', 'clinvar_mut', 'af', 'name'], 
                     sep='\t', dtype={'chrom':str})
indices = (0.999 > snp_df['af']) & (snp_df['af'] > .001)
snp_df = snp_df[indices]
snp_df.reset_index(inplace=True, drop=True)
snp_df['chrom'] = ['chr%s'%n for n in snp_df['chrom']]
snp_df.rename(columns={'chrom':'chromosome'}, inplace=True)

# Restrict to chromosomes of interest
indices = snp_df['chromosome'].isin(chromosomes)
snp_df = snp_df[indices]

# Positions need to be reduced by 1
snp_df['pos'] -= 1

# Make df to store results in
annotated_snp_df = snp_df.copy()
annotated_snp_df['pos_in_ss'] = None
annotated_snp_df['ss_ref'] = '.'
annotated_snp_df['ss_alt'] = '.'
annotated_snp_df['ss_start'] = None
annotated_snp_df['ss_stop'] = None
annotated_snp_df['ss_strand'] = None
annotated_snp_df['in_ss'] = False

# Iterate over chromosomes
for chromosome in chromosomes:
    ss_chr_df = ss_df[ss_df['chromosome'] == chromosome].copy()
    snp_chr_df = snp_df[snp_df['chromosome'] == chromosome].copy()
    
    print('chromosome: %s,\tnum ss: %d,\tnum snps: %d'%(chromosome,len(ss_chr_df),len(snp_chr_df)))

    # Sort indices
    ss_chr_df.sort_values(by='start',inplace=True)
    snp_chr_df.sort_values(by='pos',inplace=True)
    # Make snp copy with space for splice site annotation
    
    # Create ss iteratior
    ss_iterr = ss_chr_df.iterrows()
    ss_index, ss_row = next(ss_iterr)
    
    # Create SNP iterator
    snp_iterr = snp_chr_df.iterrows()
    snp_index, snp_row = next(snp_iterr)
    
    try:
        while True:
            if snp_row['pos'] < ss_row['start']:
                snp_index, snp_row = next(snp_iterr)
            elif ss_row['stop'] <= snp_row['pos']:
                ss_index, ss_row = next(ss_iterr)
            else:
  
                #print('snp in ss found at chromosome %s, pos %d'%(chromosome, snp_row['pos']))
                annotated_snp_df.loc[snp_index,'ss_start'] = ss_row['start']
                annotated_snp_df.loc[snp_index,'ss_stop'] = ss_row['stop']
                annotated_snp_df.loc[snp_index,'ss_strand'] = ss_row['strand']
                annotated_snp_df.loc[snp_index,'in_ss'] = True

                if ss_row['strand'] == '+':
                    annotated_snp_df.loc[snp_index,'pos_in_ss'] = snp_row['pos'] - ss_row['start']
                    annotated_snp_df.loc[snp_index,'ss_ref'] = snp_row['ref']
                    annotated_snp_df.loc[snp_index,'ss_alt'] = snp_row['alt']
                elif ss_row['strand'] == '-':
                    annotated_snp_df.loc[snp_index,'pos_in_ss'] = 8 - (snp_row['pos'] - ss_row['start']) 
                    annotated_snp_df.loc[snp_index,'ss_ref'] = rc(snp_row['ref'])
                    annotated_snp_df.loc[snp_index,'ss_alt'] = rc(snp_row['alt'])
                        
                # Iterate snp
                snp_index, snp_row = next(snp_iterr)
                
    except StopIteration:
        pass
                
print('Done!')
    
annotated_snp_df = annotated_snp_df[annotated_snp_df['in_ss']]
annotated_snp_df.reset_index(inplace=True,drop=True)
del annotated_snp_df['in_ss']
del annotated_snp_df['ref']
del annotated_snp_df['alt']

# Load human fasta
hg = Fasta(hg_fasta_file)

# Allocate memory
print('Allocating memory...')
annotated_snp_df['ss'] = 'N'*9

# Load all splice site sequences
print('Looking up splice site sequences...')
start_time = time.time()
chromosome_prev = None
for index, row in annotated_snp_df.iterrows():
    start = row['ss_start']
    end = row['ss_stop']
    chromosome = row['chromosome']
    if chromosome_prev != chromosome:
        chromosome_prev = chromosome
        chrom = hg[chromosome]    
        print('Working on chromosome %s...' % chromosome)
    seq = str(chrom[start:end].upper())
    if row['ss_strand']=='-':
        seq = rc(seq)
    annotated_snp_df.loc[index,'ss'] = seq
    
# Record whether splice splice site has G[TC] in right location
ss_pattern = re.compile('...G[TC]....')
indices = np.array([bool(re.match(ss_pattern,seq)) for seq in annotated_snp_df['ss']])
print('Removing %s of %s splice sites that do not have a GY at +1-2'%
      (sum(~indices), len(annotated_snp_df)))
annotated_snp_df = annotated_snp_df[indices]
 

annotated_snp_df['match'] = [r['ss_ref']==r['ss'][r['pos_in_ss']] for i,r in annotated_snp_df.iterrows()]
annotated_snp_df.head()
sum(~annotated_snp_df['match'])
num_matches = sum(annotated_snp_df['match'])
num_snps = len(annotated_snp_df)
print('%d of %d SNP refernce nucleotides match positions in splice sites.'%\
      (sum(annotated_snp_df['match']),len(annotated_snp_df)))

# Save annotated splice sites
annotated_snp_df.to_csv(exac_snp_file,index=False,sep='\t')
