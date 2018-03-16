#!/usr/bin/env python
import os, sys
import pandas as pd
import numpy as np
import pdb

from code import metadata
from code import utils

metadata_file = sys.argv[1]
group = sys.argv[2]
use_sample_data = eval(sys.argv[3])
experiment = sys.argv[4]

# Load metadata
metadata.load(metadata_file, group=group, use_sample_data=use_sample_data)

# Get files corresponding to experiment
df = metadata.get_sheet('experiments')
df.set_index('experiment',inplace=True)

# Make sure experiment is in library
assert experiment in df.index.values,\
    'Error: experiment %s not listed in metadata.'%experiment

library = df['library'][experiment]
exon_sample = df['exon_bc_sample'][experiment]
total_sample = df['total_bc_sample'][experiment]

cipher_elements = {'library':library}
library_file = metadata.encode_filename(cipher_elements,'cipher')
assert os.path.isfile(library_file),\
    'Error: file %s does not exist.'%library_file

exon_elements = {'sample':exon_sample}
exon_file = metadata.encode_filename(exon_elements,'counts')
assert os.path.isfile(exon_file),\
    'Error: file %s does not exist.'%exon_file

total_elements = {'sample':total_sample}
total_file = metadata.encode_filename(total_elements,'counts')
assert os.path.isfile(total_file),\
    'Error: file %s does not exist.'%total_file

# Validate output file
out_elements = {'experiment':experiment}
out_file = metadata.encode_filename(out_elements,'results')

# Create name of report file
report_file = metadata.filename_to_report_filename(out_file)

# Load library data frame
library_df = pd.read_csv(library_file,delim_whitespace=True,header=0)
library_df = library_df[['bc','ss','ct','otherct']].\
                rename(columns={'ct':'lib_ct','otherct':'mis_ct'})

# Load exon inclusion data frame
exon_df = pd.read_csv(exon_file,delim_whitespace=True,header=0)
exon_df = exon_df[['bc','ct']].rename(columns={'ct':'ex_ct'})

# Load total rna data frame
total_df = pd.read_csv(total_file,delim_whitespace=True,header=0)
total_df = total_df[['bc','ct']].rename(columns={'ct':'tot_ct'})

# Gather into results_df. Sort on total counts
results_df = library_df.\
            merge(exon_df,on='bc',how='left').\
            merge(total_df,on='bc',how='left')
results_df=results_df[['tot_ct','ex_ct','lib_ct','mis_ct','ss','bc']]
results_df.fillna(0,inplace=True)
results_df.sort('tot_ct',inplace=True,ascending=False)
results_df.reset_index(inplace=True,drop=True)

# Write results_df to output
results_df.to_csv(out_file,sep='\t',float_format='%d', na_rep='NaN')

# Should I write a report file? 
num_bc = len(results_df['bc'].values)
num_unique_bc = len(set(results_df['bc'].values))
num_ss = len(set(results_df['ss'].values))
assert ('9nt' in out_file) or ('11nt' in out_file),\
    'Error: cant determine the splice site length'
num_possible_ss = 2*4**7 if '9nt' in out_file else 2*4**9 
pct_ss = 100*num_ss/num_possible_ss

mean_lib_ct = results_df['lib_ct'].mean()

# Compile report
out_string =  'library_file:\t%s\n'%library_file
out_string += 'total_file:\t%s\n'%total_file
out_string += 'exon_file:\t%s\n'%exon_file
out_string += 'out_file:\t%s\n'%out_file
out_string += 'num bc:\t%d (%d unique)\n'%(num_bc,num_unique_bc)
out_string += 'num ss:\t%d (%.1f%% coverage)\n'%(num_ss,pct_ss)
for col in ['tot_ct','ex_ct','lib_ct','mis_ct']:
    out_string += '%s:\t%.1f += %.1f\n'%\
        (col, results_df[col].mean(), results_df[col].std())
out_string += '-------------\n'

# Write report
with open(report_file,'w') as f:
    f.write(out_string)