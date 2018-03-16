#!/usr/bin/env python
from __future__ import division
import os, sys
import re
import pdb
import pandas as pd

from code import metadata
from code import utils

metadata_file = sys.argv[1]
group = sys.argv[2]
use_sample_data = eval(sys.argv[3])
LID = int(sys.argv[4])
split = sys.argv[5]

# Load metadata
metadata.load(metadata_file, group=group, use_sample_data=use_sample_data)

# Get name of input file
in_filetype = 'reads'
in_elements = {'LID':LID,'split':split}
in_file = metadata.encode_filename(in_elements, in_filetype)

# Check that input files exist
assert os.path.isfile(in_file),\
    'Error: could not find input file %s.'%in_file

# Get name of input file
out_filetype = 'features'
out_elements = {'LID':LID,'split':split}
out_file = metadata.encode_filename(out_elements, out_filetype)

# Create temporary file to write to line-by-line
tmp_file = out_file+'.tmp'
tmp_file_handle = open(tmp_file,'w')

# Get name of report file
report_file = metadata.filename_to_report_filename(out_file)

# Class to track parsing failures
fail = utils.FailureTracker()

# Parse input dataframe row-by-row; fill output dataframe
in_df = pd.read_csv(in_file,delim_whitespace=True,header=0)
features_dict = {'bc':None,'ss':None,'jct':None}
#tmp_df = pd.DataFrame(columns=['sample','bc','ss','jct'])

columns=['sample','bc','ss','jct']
tmp_file_handle.write('\t'.join(columns)+'\n')

valid_samples = metadata.get_sheet('samples')['sample'].values
for i, row in in_df.iterrows():
    sample = row['sample']
    read1_sequence = row['read1_cropped']
    read2_sequence = row['read2_cropped']
    f_dict = features_dict.copy()

    if fail.test(not sample in valid_samples,'sample not vaid.'): continue;

    read1_features = metadata.parse_read_given_sample(sample,1,read1_sequence)
    if fail.test(not read1_features,'read1 not parsed.'): continue;

    read2_features = metadata.parse_read_given_sample(sample,2,read2_sequence)
    if fail.test(not read2_features,'read2 not parsed.'): continue;

    f_dict['sample'] = sample
    for key in features_dict:
        f1 = read1_features[key]
        f2 = read2_features[key]
        if fail.test((f1 and f2) and (f1!=f2),'%s differs between reads'%key):
            continue
        elif f1:
            f_dict[key] = f1
        elif f2:
            f_dict[key] = f2
        else:
            f_dict[key] = None

    # Add row to output dataframe
    #tmp_df.loc[tmp_df.shape[0]] = pd.Series(f_dict)

    # Write results line by line
    line = '\t'.join([str(f_dict[c]) for c in columns])+'\n'
    tmp_file_handle.write(line)

# Read in outfile
tmp_file_handle.close()

# Read in outfile
tmp_df = pd.read_csv(tmp_file,delim_whitespace=True)

# Count identical rows 
cols = list(tmp_df.columns.values)
tmp_df['ct'] = 1
out_df = tmp_df.groupby(cols, as_index=False).sum()

# Reorder columns and sort by counts
out_df = out_df[['ct','sample','ss','bc','jct']]
out_df.sort('ct',inplace=True,ascending=False)
out_df.reset_index(inplace=True,drop=True)

# Write output dataframe
out_df.to_csv(out_file,sep='\t')

# Compile report
num_reads_in = len(in_df)
num_reads_out = out_df['ct'].sum()
out_string = 'in_file:  %s\n'%in_file
out_string += 'out_file: %s\n'%out_file
out_string += 'num_reads_in:  %d\n'%num_reads_in
out_string += 'num_reads_out: %d\n'%num_reads_out
out_string += '%02.2f%%: Success.\n'%(100.*num_reads_out/num_reads_in)
fail_dict = fail.get_results()
for key in fail_dict:
    num_fail_type = fail_dict[key]
    out_string += '%02.2f%%: Failure, %s\n'%\
        (100.*num_fail_type/num_reads_in,key)
out_string += '-------------\n'

# Write report
with open(report_file,'w') as report_file_handle:
    report_file_handle.write(out_string)





