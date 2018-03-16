#!/usr/bin/env python

import pandas as pd
import string
import sys
from os.path import isfile
from datetime import datetime
import re
import pdb

sys.path.append('../')
from code import metadata

# Get user input
metadata_file = sys.argv[1]
group = sys.argv[2]
use_sample_data = eval(sys.argv[3])

# Print date and time
print 'Running validate_features.py on %s'%str(datetime.now())

# Function to compute the reverse complement of a sequence
complement = string.maketrans(u'ATCGN', u'TAGCN')
def rc(seq):
    seq = seq.encode('ascii','ignore') # Needed to work w/ unicod strings
    return seq.upper().translate(complement)[::-1]

# Test rc function
seq = u'ACGTNacgtn'
seq_rc = 'NACGTNACGT'
assert(rc(seq)==seq_rc)

# Verify two arguments provided
num_arguments = len(sys.argv)-1
assert num_arguments==3, 'Error: validate_features.py takes 2 argument (the metadata file name). %d arguments provided.'%num_arguments

# Verify metadata file
assert isfile(metadata_file), 'Error: metadata_file %s does not exist.'%metadata_file

# Load metadata
metadata.load(metadata_file, group=group, use_sample_data=use_sample_data)

# Extract all sheets as separate dataframes
runs_df = metadata.get_sheet('illumina_runs')
samples_df = metadata.get_sheet('samples')
experiments_df = metadata.get_sheet('experiments')
libraries_df = metadata.get_sheet('libraries')
bc_samples_df = metadata.get_sheet('bc_samples')
jct_samples_df = metadata.get_sheet('jct_samples')
amplicons_df = metadata.get_sheet('amplicons')
reads_df = metadata.get_sheet('reads')

#
# Validate shared columns
#
print '''
############################
Validating shared columns
'''

print 'Validating LIDs listed in samples sheet...'
assert set(runs_df['LID']) >= set(samples_df['LID'])

print 'Validating amplicons listed in reads sheet...'
assert set(amplicons_df['amplicon']) >= set(reads_df['amplicon'])

print 'Validating amplicons listed in samples sheet...'
assert set(amplicons_df['amplicon']) >= set(samples_df['amplicon'])

print 'Validating formats listed in amplicons sheet...'
valid_formats = [u'ssbc',u'bc',u'jct']
assert set(valid_formats) >= set(amplicons_df['format'])

print 'Validating read_nums listed in reads sheet...'
valid_read_nums = [1,2]
assert set(valid_read_nums) >= set(reads_df['read_num'])

print 'OK'

#
# Validate illumina_run files
#
print '''
############################
Validating illumina_run files
'''

for i, row in runs_df.iterrows():
    lid = row['LID']
    lid_str = '%d'%lid
    directory = metadata.get_filetype_directory('illumina_run')
    read1_file = directory + '/' + row['read1_file']
    read2_file = directory + '/' + row['read2_file']

    print 'Checking LID %d read1 file: %s...'%(lid,read1_file)
    assert lid_str in read1_file
    assert 'R1' in read1_file
    assert 'fastq.gz' in read1_file
    assert isfile(read1_file)

    print 'Checking LID %d read2 file: %s...'%(lid,read2_file)
    assert lid_str in read2_file
    assert 'R2' in read2_file
    assert 'fastq.gz' in read1_file
    assert isfile(read2_file)

print 'OK'


#
# Validate sequence reads and regular expressions
#
print '''
############################
Validating sequence reads'''

for i, row in reads_df.iterrows():
    sequence = str(row['sequence'])
    read_num = row['read_num']
    amplicon = row['amplicon']
    read_name = row['read']
    features = re.split('\W+',row['features'])

    # Check read name
    assert read_name == '%s_read%d'%(amplicon,read_num),\
        'Error: invalid read name %s.'%read_name

    # Check that read sequence parses
    parsed_dict = metadata.parse_read(read_name,sequence)
    assert parsed_dict,\
        'Error: could not parse sequence for read %s'%read_name

    # Check that all features are parsed
    for feature in features:
        assert feature in parsed_dict,\
            'Error: feature %s not found in parsed_dict %s'%\
            (feature, str(parsed_dict))
        assert len(parsed_dict[feature]) > 0,\
            'Error: feature %s is empty; should be parsed from read %s'%\
            (feature,read_name)

print 'OK'
