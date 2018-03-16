#!/usr/bin/env python
import sys
import os

print('Running compute_rnahybrid.py ...')

sys.path.append('src')
import utils

u1_fasta = 'data/U1_snRNA.fa'
all_ss_fasta = 'data/ss_GT_9nt.fa'
out_predictions = 'output/rnahybrid_predictions.txt'

# Create all 9nt GT splice sites and save to fasta
iupac_motif = 'NNNGTNNNN'
dna_seqs = utils.make_iupac_seqs(iupac_motif)
rna_seqs = [s.replace('T', 'U') for s in dna_seqs]
with open(all_ss_fasta, 'w') as f:
    for n, s in enumerate(rna_seqs):
        f.write('> %s\n%s\n' % (s, s))

# Run RNA hybrid on all GT splice sites; grab only ss and energy output
command = ('RNAhybrid -t %s -q %s -s 3utr_fly -c | ' +
           'awk \'BEGIN{FS=":"} {print $1 "\\t" $5}\' > %s') %\
          (all_ss_fasta, u1_fasta, out_predictions)
os.system(command)

print('Done!')