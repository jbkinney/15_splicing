from __future__ import division
import os
import shutil
import pandas as pd
import numpy as np
import utils
import re
import pdb
import utils
import time
import sys
import glob

import metadata

def clean_results_and_samples():
    ''' WARNING: HARD-CODED DIRECTORIES HERE '''

    # Clean results directory
    d = 'results'
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.mkdir(d)

    # Clean sampled illumina files
    files = glob.glob('data/illumina_runs/sample*.gz')
    for file in files:
        os.remove(file)

def clean_intermediate_and_tmp_directories():
    # Get information about directory structure and which dirs to clean
    filetypes_df = metadata.get_sheet('pipeline_filetypes')
    dir_col = filetypes_df['directory']
    clean_col = filetypes_df['clean_at_start']

    # Check the existence of directories that won't be cleaned
    dirs_to_check = dir_col[clean_col==False].values
    for d in dirs_to_check:
        assert os.path.isdir(d),\
        'Error: directory %s doesnt exist.'%d

    # Make a list of directories meant for cleaning
    dirs_to_clean = dir_col[clean_col==True].values
    for d in dirs_to_clean:

        # If direcotry exists, remove it
        if os.path.isdir(d):
            shutil.rmtree(d)

        # Make directory
        #print '\t Creating directory %s'%d
        os.mkdir(d)

    # Clean tmp dir
    utils.clean_dir(utils.TMP_DIR)

def validate_metadata():
    '''
    Makes sure that the metadata file satisfies all consistency requirements.
    '''
    metadata_file, group = metadata.get_filename_and_group()
    use_sample_data = metadata.get_use_sample_data()
    out_file = 'data/report.metadata.txt'
    command = './scripts/validate_metadata.py %s %s %s > %s\n'%\
        (metadata_file, group, use_sample_data, out_file)

    # Dont use cluster for this
    utils.submit_and_complete_jobs([command],
        'validate_metadata',use_cluster=False)


def validate_filename_handling():
    '''
    Makes sure that filename handling functions in metadata pass 
    a series of tests.
    '''
    metadata_file, group = metadata.get_filename_and_group()
    tests_file = 'data/tests.xlsx'
    out_file = 'data/report.tests.txt'
    command = './scripts/validate_filename_handling.py %s %s %s > %s\n'%\
        (metadata_file, group, tests_file, out_file)

    # Dont use cluster for this
    utils.submit_and_complete_jobs([command],
        'validate_filename_handling',use_cluster=False)


def validate_illumina_runs():
    '''
    Makes sure that all illumina run files listed in metadata exist. 
    '''
    directory = metadata.get_filetype_directory('illumina_run')
    df = metadata.get_sheet('illumina_runs')
    for i, row in df.iterrows():

        # Validate read1 file
        read1_file_name = directory + '/' + row['read1_file']
        #print '\t Verifying %s'%read1_file_name
        assert os.path.isfile(read1_file_name),\
            'Error: could not find illumina_run file %s.'%read1_file_name

        # Validate read2 file
        read2_file_name = directory + '/' + row['read2_file']
        #print '\t Verifying %s'%read2_file_name
        assert os.path.isfile(read2_file_name),\
            'Error: could not find illumina_run file %s.'%read2_file_name


def sample_illumina_run_files(num_sampled_reads):
    '''
    Samples illumina_run files. Use to troubleshoot pipeline.
    NOTE: Run this *before* validate_illumina_runs
    '''
    num_sampled_lines = 4*num_sampled_reads
    directory = metadata.get_filetype_directory('illumina_run')

    # Remove existing sample files
    df = metadata.get_sheet('pipeline_filetypes')
    df.set_index('file_type',inplace=True)
    base_glob = df['glob']['illumina_run']
    files_glob = directory + '/' + base_glob
    files_to_remove = glob.glob(files_glob)
    for f in files_to_remove:
        assert 'sample.' in f,\
            'Error: file %s does not appear to be a sampled illumina_run.'%f
        os.remove(f)

    df = metadata.get_sheet('illumina_runs')
    basenames = list(df['read1_file']) + list(df['read2_file'])
    out_files = [directory + '/' + b for b in basenames]
    scripts = []
    for out_file in out_files:

        p = re.compile('.*/sample\.(.*)')
        m = re.match(p,out_file)
        in_file = directory + '/' + m.group(1)
        try:
            assert os.path.isfile(in_file),\
                'Error: file %s does not exist.'%in_file
        except:
            pdb.set_trace()

        command = 'gzip -cd %s | head -n %d | gzip > %s\n'%\
            (in_file, num_sampled_lines, out_file)
        scripts.append(command)

    # Use cluster
    utils.submit_and_complete_jobs(scripts,
        'sample_illumina_run_files',use_cluster=True)


def make_split_fastq_files(reads_per_split):
    '''
    Splits an illumina_run file, specified by LID and read_num, into
    multiple split_fastq files.

    input file_type: illumina_runs
    output file_type: split_fastq 
    '''
    in_file_type = 'illumina_run'
    
    # Clean output directory
    out_file_type = 'split_fastq'
    out_dir = metadata.get_filetype_directory(out_file_type)
    utils.clean_dir(out_dir)

    LIDs = metadata.get_all_element_values('LID')
    read_nums = metadata.get_all_element_values('read_num')

    commands = []
    for LID in LIDs:
        for read_num in read_nums:

            # Get in file
            try:
                files = metadata.get_existing_filenames(in_file_type,\
                    {'LID':LID,'read_num':read_num})
            except:
                pdb.set_trace()
            assert len(files)>0,\
                'Error: no %s files matching criteria.'%in_file_type
            assert len(files)==1,\
                'Error: finding multiple %s files matching criteria.'%in_file_type
            in_file_name = files[0]

            num_lines = int(4*reads_per_split)
            d = {'LID':LID,'read_num':read_num,'split':''} 
            out_file_prefix = metadata.encode_filename(d,out_file_type)
            command = 'zcat %s | split -a 3 -l %d - %s'%\
                (in_file_name,num_lines,out_file_prefix)
            commands.append(command)

    # Execute commands; use cluster
    utils.submit_and_complete_jobs(commands,
        'make_split_fastq_files',use_cluster=True)
  

def make_reads_files(): 
    '''
    Given a read1_split_fastq and read2_fastgz files, specified by LID and split, 
    this function creates a reads file and corresponding report. 
    '''

    # Clean output directory
    out_file_type = 'reads'
    out_dir = metadata.get_filetype_directory(out_file_type)
    utils.clean_dir(out_dir)

    metadata_file, group = metadata.get_filename_and_group()
    use_sample_data = metadata.get_use_sample_data()
    LIDs = metadata.get_all_element_values('LID')
    commands = []
    for LID in LIDs:
        splits = metadata.get_splits_from_LID(LID)
        for split in splits:
            command = './scripts/make_reads_files.py %s %s %s %s %s'%\
                (metadata_file, group, use_sample_data, LID, split)
            commands.append(command)

    # Execute commands; use cluster for this
    utils.submit_and_complete_jobs(commands,
        'make_reads_files',use_cluster=True)


def make_features_files(): 
    '''
    For each reads file , this function creates
    a features file and a corresponding report.
    '''
    
    # Clean output directory
    out_file_type = 'features'
    out_dir = metadata.get_filetype_directory(out_file_type)
    utils.clean_dir(out_dir)

    metadata_file, group = metadata.get_filename_and_group()
    use_sample_data = metadata.get_use_sample_data()

    LIDs = metadata.get_all_element_values('LID')
    commands = []
    for LID in LIDs:
        splits = metadata.get_splits_from_LID(LID)
        for split in splits:
            command = './scripts/make_features_files.py %s %s %s %s %s'%\
                (metadata_file, group, use_sample_data, LID, split)
            commands.append(command)

    # Execute commands; use cluster for this
    utils.submit_and_complete_jobs(commands,
        'make_features_files',use_cluster=True)


def make_efficiency_file():
    '''
    Processes report files assess pipeline performance
    '''

    metadata_file, group = metadata.get_filename_and_group()

    out_file = 'results/efficiency.%s.txt'%group

    report_pattern = re.compile(\
        'num_reads_in:\s+(?P<reads_in>[0-9]+)\s+'+ \
        'num_reads_out:\s+(?P<reads_out>[0-9]+)\s+.*')

    columns = ['LID','split','raw_reads','sorted_reads','parsed_reads']
    tmp_df = pd.DataFrame(columns=columns)

    # Get name of all reads files files
    reads_files = metadata.get_existing_filenames('reads')
    for n, reads_file in enumerate(reads_files):
        elements_dict = metadata.decode_filename(reads_file,'reads')

        # Parse report file
        reads_report_file = metadata.encode_filename(\
            elements_dict,'reads',report=True)
        assert os.path.isfile(reads_report_file), \
            'Could not find reads report file %s'%reads_report_file
        with open(reads_report_file,'r') as f:
            content = f.read()
        m = re.search(report_pattern,content)
        try:
            assert m, 'Could not parse report file %s.'%reads_report_file
        except:
            pdb.set_trace()
        reads_dict = m.groupdict()

        # Parse features file
        features_report_file = metadata.encode_filename(\
            elements_dict,'features',report=True)
        assert os.path.isfile(features_report_file), \
            'Could not find features report file %s'%features_report_file
        with open(features_report_file,'r') as f:
            content = f.read()
        m = re.search(report_pattern,content)
        assert m, 'Could not parse report file %s.'%features_report_file
        features_dict = m.groupdict()

        # Record results
        assert int(reads_dict['reads_out'])==int(features_dict['reads_in']),\
            'Number of reads out from reads file does not match'+\
            'reads in from features file.'

        line_dict = {
            'LID':elements_dict['LID'],
            'split':elements_dict['split'],
            'raw_reads':int(reads_dict['reads_in']),
            'sorted_reads':int(reads_dict['reads_out']),
            'parsed_reads':int(features_dict['reads_out'])
        }
        tmp_df.loc[n,:] = line_dict

    # Marginalize over splits
    out_df = tmp_df.groupby('LID').sum()
    del out_df['split']

    # Comptue efficiencies
    out_df['frac_sorted'] = out_df['sorted_reads']/out_df['raw_reads']
    out_df['frac_parsed'] = out_df['parsed_reads']/out_df['raw_reads']

    # Annotate LIDs with group and description
    df = metadata.get_sheet('illumina_runs')
    df.set_index('LID',inplace=True,drop=True)
    out_df['group'] = group
    out_df['description'] = df['description']

    # Write data frame to distk
    out_df.to_csv(out_file,sep='\t')


def make_counts_files(): 
    '''
    For each sample, this function collates counts for all sets
    of features across all relevant features files (specified by corresponding
    LID).
    '''

    out_file_type = 'counts'
    out_dir = metadata.get_filetype_directory(out_file_type)
    utils.clean_dir(out_dir)

    metadata_file, group = metadata.get_filename_and_group()
    use_sample_data = metadata.get_use_sample_data()

    samples = metadata.get_all_element_values('sample')
    commands = []
    for sample in samples:
        command = './scripts/make_counts_files.py %s %s %s %s'%\
            (metadata_file, group, use_sample_data, sample)
        commands.append(command)

    # Execute commands, use cluster for this
    utils.submit_and_complete_jobs(commands,
        'make_counts_files',use_cluster=True)


def make_cipher_files(): 
    '''
    For each library, output a cipher that maps splice sites to barcodes. 
    '''

    out_file_type = 'cipher'
    out_dir = metadata.get_filetype_directory(out_file_type)
    utils.clean_dir(out_dir)

    metadata_file, group = metadata.get_filename_and_group()
    use_sample_data = metadata.get_use_sample_data()

    libraries = metadata.get_all_element_values('library')
    commands = []
    for library in libraries:
        command = './scripts/make_cipher_files.py %s %s %s %s'%\
            (metadata_file, group, use_sample_data, library)
        commands.append(command)

    # Execute commands, use cluster for this
    utils.submit_and_complete_jobs(commands,
        'make_cipher_files',use_cluster=True)


def make_results_files():
    '''
    Given an experiment, collect counts for the exon_bc_sample and 
    total_bc_sample, as well as information about the library.
    Do NOT clean out directory before running
    WARNING: REFERENCES HARD-CODED DIRECTORIES
    '''

    out_file_type = 'results'
    out_dir = metadata.get_filetype_directory(out_file_type)

    metadata_file, group = metadata.get_filename_and_group()
    use_sample_data = metadata.get_use_sample_data()

    experiments = metadata.get_all_element_values('experiment')
    commands = []
    for experiment in experiments:
        command = './scripts/make_results_files.py %s %s %s %s'%\
            (metadata_file, group, use_sample_data, experiment)
        commands.append(command)

    # Execute commands, use cluster for this
    utils.submit_and_complete_jobs(commands,
        'make_results_files',use_cluster=True)

    # Finally, copy ciphers and junction counts to results folder
    os.system('cp intermediate/ciphers/* results/.')
    os.system('cp intermediate/counts/*_jct_* results/.')
    