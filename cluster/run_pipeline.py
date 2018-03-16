#!/usr/bin/env python

from __future__ import division
import os
import pandas as pd
import numpy as np
from glob import glob
import shutil
import pdb
import sys
import time

#from code import pipeline
from code import metadata
from code import utils
from code import pipeline

def run_stage(stage_num, text):
    '''
    Says whether or not to run stage with given stage_num.
    Prints description if so
    '''
    val = (FIRST_STAGE <= stage_num) and (LAST_STAGE >= stage_num)
    if val:
        print '\n##### Stage %d: %s... #####'%(stage_num,text)
    return val

###############################################################################
### USER OPTIONS ###############################################################################

# Flag for whether or not to use sampled data
metadata_file = 'data/metadata.xlsx'
groups_to_run = ['brca2_9nt','ikbkap_9nt','smn1_9nt','brca2_11nt']
start_over = True
use_sample_data = False
num_sampled_reads = int(1E6)
reads_per_split = int(1E6)

FIRST_STAGE = -1
LAST_STAGE = 10
USE_CLUSTER = True

###############################################################################
###############################################################################

# Set start time
start_time = time.time()

# Whether or not to use cluster and, if so, check for jobs currently running
utils.use_cluster(USE_CLUSTER)
if USE_CLUSTER:
    # Check for existing jobs
    qstat_file = '%s/qstat_output.txt'%utils.TMP_DIR
    utils.execute_at_commandline('qstat > %s'%qstat_file)
    assert os.path.isfile(qstat_file)
    assert os.path.getsize(qstat_file)==0,\
        'Error: jobs are still running on the cluster!'
    os.remove(qstat_file)

# If start over, clean results and sampled files
if start_over and run_stage(-1,'Cleaning results and samples.'):
    pipeline.clean_results_and_samples()

# Run pipeline over all specified groups
for group in groups_to_run:

    # Set start time
    group_start_time = time.time()

    # Load metadata for specified group
    print 'Processing group %s using %s...'%(group,metadata_file)
    metadata.load(metadata_file, group=group, use_sample_data=use_sample_data)

    # Clean directories / make directory structure
    if run_stage(0,'Cleaning pipeline directories for group %s'%group):
        pipeline.clean_intermediate_and_tmp_directories()

    # Create sample data if instructed to
    if use_sample_data:
        if run_stage(1,'Sampling illumina_run_files for group %s'%group):
            pipeline.sample_illumina_run_files(num_sampled_reads)

    # Vaidate metdata files 
    # Note: checks for existence of sampled data if sampling, so musht come
    # after satge 1
    if run_stage(2,'Validating metadata file %s for group %s'%\
        (metadata_file,group)):
        pipeline.validate_metadata()

    # Validate illumina run files
    if run_stage(3,'Validating illumina_run files for group %s'%group):
        pipeline.validate_illumina_runs()

    # Map illumina_runs -> split_fastq
    if run_stage(4,'Creating split_fastq files for group %s'%group):
        pipeline.make_split_fastq_files(reads_per_split)

    # Map split_fastq -> reads
    if run_stage(5,'Creating reads files for group %s'%group):
        pipeline.make_reads_files()

    # Map reads -> features
    if run_stage(6,'Creating features files for group %s'%group):
        pipeline.make_features_files()

    # Determine pipeline efficiency
    if run_stage(7,'Creating efficiency file for group %s'%group):
        pipeline.make_efficiency_file()

    # Map features -> counts
    if run_stage(8,'Creating counts files for group %s'%group):
        pipeline.make_counts_files()

    # Map counts -> cipher
    if run_stage(9,'Creating cipher files for group %s'%group):
        pipeline.make_cipher_files()

    # Map counts, ciphers -> results
    if run_stage(10,'Creating results files for group %s'%group):
        pipeline.make_results_files()

    # Announce group run time
    group_end_time = time.time()
    group_run_time_min = (group_end_time-group_start_time)/60.0
    print 'Group %s done! Took %.1f minutes.'%\
        (group,group_run_time_min)

# Announce full runtime
end_time = time.time()
run_time_min = (end_time-start_time)/60.0
print 'Entire pipeline took %.1f minutes. Groups processed: %s'%\
    (run_time_min,', '.join(groups_to_run))



