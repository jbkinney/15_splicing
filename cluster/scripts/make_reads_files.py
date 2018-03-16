#!/usr/bin/env python
import os, sys
import pandas as pd
import time
import pdb

from code import metadata
from code import utils

# Get user input
metadata_file = sys.argv[1]
group = sys.argv[2]
use_sample_data = eval(sys.argv[3])
LID = sys.argv[4]
split = sys.argv[5]

# Start timer
total_time_start = time.time()

# Load metadata
metadata.load(metadata_file, group=group, use_sample_data=use_sample_data)

# Get inpt file names
in_file_type = 'split_fastq'
read1_file = metadata.encode_filename(
    {'LID':LID, 'split':split, 'read_num':1},in_file_type)
read2_file = metadata.encode_filename(
    {'LID':LID, 'split':split, 'read_num':2},in_file_type)

# Check that input files exist
assert os.path.isfile(read1_file),\
    'Error: could not find input file %s.'%read1_file
assert os.path.isfile(read2_file),\
    'Error: could not find input file %s.'%read2_file

# Get name of output file
out_file_type = 'reads'
out_elements = {'LID':LID,'split':split}
out_file = metadata.encode_filename(out_elements, out_file_type)
out_file_handle = open(out_file,'w')

# Get name of report file
report_file = metadata.filename_to_report_filename(out_file)
report_file_handle = open(report_file,'w')

##################
# Below is adapted from parse_reads_9nt.py

# Get LID
elements1 = metadata.decode_filename(read1_file,in_file_type)
elements2 = metadata.decode_filename(read2_file,in_file_type)
assert elements1['LID'] == elements2['LID']
LID = elements1['LID']

# Get barcodes and corresponding samples. Store as dict.
barcode_to_sample = {}
df = metadata.get_sheet('samples')
LID_df = df[df['LID']==LID]
for i,row in LID_df.iterrows():
    sample = row['sample']
    barcode = row['barcode']
    barcode_to_sample[barcode] = sample

# Class to track parsing failures
fail = utils.FailureTracker()

# Process sequences one-by-one: This is where the main processing happens
r1_f = open(read1_file)
r2_f = open(read2_file)
stop = False
num_reads = 0
num_successes = 0
columns=['sample','read1_cropped','read2_cropped']
out_file_handle.write('\t'.join(columns)+'\n')
while not stop:

    # Provide feedback
    if num_reads%10000 == 0 and num_reads > 0:
        total_time = time.time() - total_time_start  
        successes_pct = 100*num_successes/num_reads

    # Get reads; halt loop if reads have length 0
    read1 = utils.get_next_read_from_fastq(r1_f)
    read2 = utils.get_next_read_from_fastq(r2_f)
    if len(read1) == 0 or len(read2) == 0:
        stop = True
        continue
    num_reads += 1

    # Determine sequence sample by matching barcode, and clip barcode
    sample1, tag_length1 = utils.match_barcode(read1,barcode_to_sample, 
        search_area=10)
    sample2, tag_length2 = utils.match_barcode(read2, barcode_to_sample, 
        search_area=10)

    # Check that both reads point to specific samples
    if fail.test(not sample1,'sample1 is not vaid.'):
        continue
    if fail.test(not sample2,'sample2 is not vaid.'):
        continue

    #Check that both samples are the same
    if fail.test(not sample1==sample2,'sample1 do not match sample2.'):
       continue

    sample = sample1 

    # Check that both reads are the right length
    assert len(read1) > tag_length1, 'Read 1 is not long enough'
    assert len(read2) > tag_length2, 'Read 2 is not long enough'
    read1_clipped = read1[tag_length1:]
    read2_clipped = read2[tag_length2:]

    out_file_handle.write('\t'.join([sample,read1_clipped,read2_clipped])+'\n')
    num_successes += 1


# Compile report
num_reads_in = num_reads
num_reads_out = num_successes
out_string = 'in_file_read1:  %s\n'%read1_file
out_string = 'in_file_read2:  %s\n'%read2_file
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

# Close report file
report_file_handle.close()