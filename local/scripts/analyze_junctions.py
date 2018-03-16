#!/usr/bin/env python
from __future__ import division
import pandas as pd
import numpy as np

import os
import sys
import re
import pdb
import glob

print('Running analyze_junctions.py ...')

sys.path.append('src')
import utils

# Specify input and output files
in_files_glob = 'from_pipeline/counts.*_jct*.txt'
out_file='output/junction_displacements.txt'

# Load junction data
jct_files = glob.glob(in_files_glob)
pattern = re.compile('counts\.(?P<name>.+)\.txt')
df_dict = {}
for f in jct_files:
    m = re.search(pattern,f)
    d = m.groupdict()
    name = d['name']
    if 'brca2' in name:
        name += '_rep1'
    experiment = utils.splice(name,[0,1,3,4])
    df = pd.read_csv(f,sep='\t',usecols=[1,4])
    df['length'] = [len(seq) for seq in df['jct']]
    df_dict[experiment] = df.copy()

# Compute splicing location (0 is expected) for each locus
expected_length_dict = {
        'brca2':31, # correct, 17.07.20
        'smn1':28, # correct, 17.07.20
        'ikbkap':31 # correct, 17.07.20
}
min_displacement = -5
max_displacement = +5
displacement_values=list(range(min_displacement,max_displacement+1,1))

# Compute displacements for each locus
experiments = list(df_dict.keys())
experiments.sort()
cols = ['<%d'%min_displacement] +\
       displacement_values + \
       ['>%d'%max_displacement]
displacement_df = pd.DataFrame(index=experiments,columns=cols)
displacement_df.index.name = 'experiment'
for experiment in experiments:
    print('Loading junction information on %s...'%experiment)
    df = df_dict[experiment]
    gene = utils.splice(experiment,[0])
    displacements = df['length']-expected_length_dict[gene]
    for d in displacement_values:
        displacement_df.loc[experiment,d]=\
            ((displacements==d)*df['ct']).sum()
    displacement_df.loc[experiment,cols[0]] = \
        sum(displacements < min_displacement)
    displacement_df.loc[experiment,cols[-1]] = \
        sum(displacements > max_displacement)
displacement_df.to_csv(out_file,sep='\t')
print('Output written to %s.'%out_file)

print('Done!')
