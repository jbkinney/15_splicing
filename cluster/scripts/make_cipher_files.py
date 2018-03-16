#!/usr/bin/env python
import os, sys
import pandas as pd
import pdb
import numpy as np

from code import metadata
from code import utils

metadata_file = sys.argv[1]
group = sys.argv[2]
use_sample_data = eval(sys.argv[3])
library = sys.argv[4]

# Load metadata
metadata.load(metadata_file, group=group, use_sample_data=use_sample_data)

# Get name of counts file
in_elements = {'sample':library}
in_file = metadata.encode_filename(in_elements,'counts')

# Make sure input file exists
assert os.path.isfile(in_file),\
    'Error: could not find input file %s.'%in_file

# Create name of cipher file
out_elements = {'library':library}
out_file = metadata.encode_filename(out_elements,'cipher')

# Create name of report file
report_file = metadata.filename_to_report_filename(out_file)

# Create empty output file
#utils.execute_at_commandline('touch %s'%out_file)

# Create empty report file
#utils.execute_at_commandline('touch %s'%report_file)

# Load input file as data frame
in_df = pd.read_csv(in_file,delim_whitespace=True,header=0)

 # Grab the most prevalent splice site for each barcode
in_df.sort(['ss','ct'],axis=0,ascending=[True,False],inplace=True)
g = in_df.groupby('bc')
out_df = g.first()

# Add 'ssnum' column: total number of associated ss per bc
tmp_df = g.size()
out_df.loc[tmp_df.index,'numss'] = tmp_df

# Add 'totct' column: total number of associated ss per bc
tmp_df = g.agg(np.sum)
out_df.loc[tmp_df.index,'totct'] = tmp_df['ct']
out_df['otherct'] = out_df['totct']-out_df['ct']

# Sort out_file
#out_df.sort('ct',axis=0,ascending=False,inplace=True)

# Compute ok indices
# Criterion: number of counts for most prevealent ss has to be
# At least 2 counts, 
# and at least 4 times the total counts for all the rest of the ss
ok_indices = (out_df['ct'] >= 2) \
    & (out_df['ct'] >= 4*out_df['otherct'])

# Get final dataframe of associated indices
out_df = out_df.loc[ok_indices,['ss','ct','otherct','numss']]
out_df.sort('ct',axis=0,ascending=False,inplace=True)

# Write output file
out_df.to_csv(out_file,sep='\t')

# Compute average number of bc per ss
num_bcs_per_ss = out_df.groupby('ss').size()
mean_bcs_per_ss = np.mean(num_bcs_per_ss)
std_bc_per_ss = np.std(num_bcs_per_ss)

# Compile report
out_string = 'in_file: %s\n'%in_file
out_string += 'out_file: %s\n'%out_file
out_string += 'bc taged: %d\n'%out_df.shape[0]
out_string += 'ss taged: %d\n'%len(num_bcs_per_ss)
out_string += 'bc per ss: %.2f +- %.2f\n'%(mean_bcs_per_ss,std_bc_per_ss)
out_string += 'reads per association: ???\n'
out_string += 'association purity: ???\n'
out_string += '-------------\n'

# Write report
with open(report_file,'w') as report_file_handle:
    report_file_handle.write(out_string)
