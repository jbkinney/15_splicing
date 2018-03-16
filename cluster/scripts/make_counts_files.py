#!/usr/bin/env python
import os, sys
import pandas as pd
import pdb

from code import metadata
from code import utils

metadata_file = sys.argv[1]
group = sys.argv[2]
use_sample_data = eval(sys.argv[3])
sample = sys.argv[4]

# Load metadata
metadata.load(metadata_file, group=group, use_sample_data=use_sample_data)

# Get LID corresponding to sample
df = metadata.get_sheet('samples')
df.set_index('sample',inplace=True)
LID = df['LID'][sample]

# Get name of all input files
in_file_type = 'features'
in_files = metadata.get_existing_filenames(in_file_type,{'LID':LID})

# Get anem of output file
out_file_elements = {'sample':sample}
out_file = metadata.encode_filename(out_file_elements,'counts')

# Get name of report file
report_file = metadata.filename_to_report_filename(out_file)

# Add relevent rows to tmp_df file-by-file
tmp_df = pd.DataFrame()
for in_file in in_files:
    
    # Load input file as data frame
    in_df = pd.read_csv(in_file,delim_whitespace=True,header=0)

    # Make sure colums are correct
    assert set(in_df.columns) >= set(['sample','ct','bc','ss','jct']),\
        'Error: Incorrect in_df columns: %s.'%str(in_df.columns)

    # Get rows of interest
    rows_df = in_df[in_df['sample']==sample]
    del rows_df['sample']

    # Append rows to tmp_df
    tmp_df = tmp_df.append(rows_df,ignore_index=True)

# Count identical rows 
cols = ['bc','ss','jct']
out_df = tmp_df.groupby(cols, as_index=False).sum()

# Reorder columns and sort by counts
out_df = out_df[['ct','ss','bc','jct']]
out_df.sort('ct',inplace=True,ascending=False)
out_df.reset_index(inplace=True,drop=True)

# Write output dataframe
out_df.to_csv(out_file,sep='\t')

# Compile report
num_reads_processed = out_df['ct'].sum()
out_string = 'in_files: \n'
for in_file in in_files:
    out_string += '\t%s\n'%in_file
out_string += 'out_file: %s\n'%out_file
out_string += 'num_reads_out: %d\n'%num_reads_processed
out_string += '-------------\n'

# Write report
with open(report_file,'w') as report_file_handle:
    report_file_handle.write(out_string)
