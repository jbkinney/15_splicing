from __future__ import division
import os
import shutil
import pandas as pd
import numpy as np
import re
import subprocess
import signal
import time
import sys
import string
import commands # Depreciated module; should replace
import glob
import pdb

TMP_DIR = 'tmp/'

USE_CLUSTER = False

WAIT_TIME = 30

def use_cluster(use_cluster):
    '''
    Set True to run commands locally. 
    Set False to run commands on cluster
    '''
    global USE_CLUSTER
    assert type(use_cluster) == bool,\
        'Error: use_cluster must be boolean.'
    USE_CLUSTER = bool(use_cluster)

def clean_dir(directory):
    '''
    Cleans a specified directory
    '''

    # Remove directory and all contents if it exists
    if os.path.isdir(directory):
        shutil.rmtree(directory)

    # Create directory
    os.mkdir(directory)

def clean_dirs(directories):
    '''
    Cleans multiple directories
    '''

    for directory in directories:
        clean_dir(directory)

def get_substring(s,pattern,group_num):
    ''' 
    get_substring:
        Returns a substring given a string s, regular expression 
        pattern, and a group number. 
    '''
    
    p = re.compile(pattern)
    m = re.match(pattern,s)
    groups = m.groups()

    assert len(groups) > 0,\
        'Error: no match to pattern %s found in string %s.'%(pattern,s)
    assert len(groups) > group_num,\
        'Error: improper group number %d'%group_num

    return groups[group_num]

def execute_at_commandline(command):
    '''
    Executes a command at the command line. Prevents the 'Broken pipe' errors
    caused by os.system()
    '''
    return subprocess.check_output(command, 
        shell=True, 
        preexec_fn=lambda:signal.signal(signal.SIGPIPE, signal.SIG_DFL))

# Useful way to give feedback
def give_feedback(feedback):
    func_name =sys._getframe(1).f_code.co_name
    print '\nIn '+func_name+': '+feedback,
    sys.stdout.flush()

class Job:
    ''' Class describing a submitted job '''
    def __init__(self,prefix,num,content,registry):
        self.num = num
        self.prefix = prefix
        self.name = '%s_%05d.sh'%(prefix,num)
        self.script_file = TMP_DIR + self.name
        self.success_file =  TMP_DIR + 'success.' + self.name
        self.error_file = self.script_file + '.e'
        self.out_file = self.script_file + '.o'
        self.script_content = content + '\necho "Success" > %s\n'%\
            self.success_file
        self.cluster_command = 'qsub -cwd -e %s -o %s ./%s' %\
            (self.error_file, self.out_file, self.script_file)
        self.local_command = './%s'%self.script_file
        self.attempts=0

        # Record job
        registry.append(self)

        # Write script file
        with open(self.script_file, 'w') as f:
            f.write(self.script_content)
        assert os.path.isfile(self.script_file)
        execute_at_commandline('chmod +x %s'%self.script_file)

        # Register status of job
        self.status = 'prepared'

    def submit(self):
        ''' Submit job to cluster and register job number '''

        # Remove any output files that were
        for f in [self.success_file, self.out_file, self.error_file]:
            if os.path.isfile(f):
                os.remove(f)

        # Submit job and get job number
        output = execute_at_commandline(self.cluster_command)
        pattern = re.compile('Your job ([0-9]+)\s')
        match = re.search(pattern,output)
        assert match, 'Error: could not parse output %s'%output
        self.job_id = int(match.group(1))
        self.status = 'submitted'
        self.attempts += 1

    def kill(self):
        self.kill_command = 'qdel %d'%self.job_id
        execute_at_commandline(self.kill_command)

    def run_locally(self):
        execute_at_commandline(self.local_command)
        self.status = 'submitted'

    def did_error_occur(self):
        return os.path.getsize(self.error_file)>0

    def do_runtime_files_exist(self):
        return os.path.isfile(self.out_file) and\
            os.path.isfile(self.error_file)

    def was_success_recorded(self):
        return os.path.isfile(self.success_file)

job_report_pattern = re.compile('\s+(?P<num>[0-9]{8})\s+.*\s+jkinney\s+(?P<status>[a-z]+)\s+.*\s+Full jobname:\s+(?P<name>\S+)')

def inspect_jobs(job_list,prefix):

    # Make jobs_dict, indexed by job name
    job_names = [j.name for j in job_list]
    job_dict = dict(zip(job_names,job_list))

    # Get and parse list of jobs on cluster
    qstat_output = commands.getoutput('qstat -r ')
    qstat_matches = re.finditer(job_report_pattern,qstat_output)

    # Record status of all jobs on cluster
    for m in qstat_matches:
        g = m.groupdict()
        name = g['name']
        status = g['status']
        job = job_dict[name]
        job.status = 'on_cluster'
        job.status_on_cluster = status

    # Get and parse list of completed jobs
    success_files = [TMP_DIR+f for f in os.listdir(TMP_DIR) \
        if prefix in f]
    job_completions = [j for j in job_list if j.success_file in success_files]
    for job in job_completions:
        job.status = 'completed'

    # Get and parse list of jobs in which errors were detected
    error_files = glob.glob(TMP_DIR+'/*.e')
    job_failures = [j for j in job_list\
        if j.error_file in error_files and os.path.getsize(j.error_file)>0]
    for job in job_failures:
        job.status = 'failed'
        with open(job.error_file,'r') as f:
            job.error_message = f.read()

# Submits a list of scripts to be run as separate jobs
# Waits for all jobs to complete before continuing
def submit_and_complete_jobs(script_content_list,prefix='script',use_cluster=True):
    global USE_CLUSTER

    
    # Both the global USE_CLUSTER and local use_cluster variables have to
    # be true for process to go to cluster
    use_cluster = use_cluster and USE_CLUSTER

    # Time code
    start_time = time.time()

    # Create set of all jobs and record in registry
    job_list = []
    num_jobs = len(script_content_list)
    give_feedback('Preparing to run %d jobs...'%(num_jobs))
    for n, content in enumerate(script_content_list):
        if n%100==0 and n>0:
            give_feedback('Preparing job number %d/%d...'%(n,num_jobs))
        job = Job(prefix=prefix,num=n,content=content,
                registry=job_list)

    # If running locally, run jobs serially
    if not use_cluster:
        for job in job_list:

            # Execute script
            give_feedback('Executing job %s...'%job.name)
            job.run_locally()

            # Make sure that success file has been written
            assert job.was_success_recorded(),\
                'Error: success file %s not written after execution of %s'%\
                (job.success_file,job.name)

    # If using cluster
    elif use_cluster:
        
        # Submit all jobs
        give_feedback('Submitting %d jobs to cluster...'%num_jobs)
        submission_start_time = time.time()
        for n, job in enumerate(job_list):
            if n%100==0 and n>0:
                give_feedback('Submitting job number %d/%d...'%(n,num_jobs))
            job.submit()
        submission_time = time.time() - submission_start_time
        give_feedback('Job submission took %d seconds, %.3f seconds per job.'%\
            (submission_time,submission_time/num_jobs))

        # Restart jobs that fail
        give_feedback('Waiting for jobs to complete...')
        job_start_time = time.time()
        max_attempts = 3
        jobs_completed = []
        jobs_remaining = job_list
        while len(jobs_completed) < len(job_list):
            
            # Give feedback
            job_run_time = time.time() - job_start_time
            give_feedback(
                '%d/%d jobs remain after %d seconds; waiting %d seconds...'%\
                (len(jobs_remaining), len(job_list), job_run_time, WAIT_TIME))

            # Wait
            time.sleep(WAIT_TIME)

            # Inspect jobs
            inspect_jobs(job_list, prefix)

            # If any jobs failed, kill and restart
            for job in job_list:
                if job.status=='failed':
                    if job.attempts>=max_attempts:
                        give_feedback(
                            'Job %s failed after %d attempts. Quitting...'%\
                            (job.name,job.attempts))
                        raise Exception
                    else:
                        give_feedback('Job %s failed. Restarting...'%job.name)
                        try:
                            job.kill()
                        except:
                            pass
                        job.submit()

            # Count jobs remaining and completed
            jobs_remaining = [j for j in job_list \
                if j.status in ['on_cluster','submitted']] 
            jobs_completed = [j for j in job_list if j.status=='completed'] 

        # Make sure all jobs completed
        for job in job_list:
            assert job.status=='completed', \
                'Error: job %s unexpectedly not completed'%job.name
            
        job_run_time = time.time() - job_start_time
        give_feedback('All %d jobs finished in %d seconds.'%\
            (len(job_list),job_run_time))

    # Announce job completion
    total_execution_time = time.time() - start_time
    give_feedback('Done. Total execution time was %d seconds.\n'%\
        total_execution_time)


# Function to compute the reverse complement of a sequence
complement = string.maketrans('ATCGN', 'TAGCN')
def rc(seq):
    return seq.upper().translate(complement)[::-1]

# Return the next read from a fastq file
def get_next_read_from_fastq(f):
    f.readline()
    read = f.readline().strip()
    f.readline()
    f.readline()
    return read

# Finds the barcode corresponding to each sequence
def match_barcode(seq,barcodes_dict,search_area=20):
    tag_length = 0
    region = False
    for barcode in barcodes_dict.keys():
        k = seq.find(barcode,0,search_area)
        if k >= 0:
            region = barcodes_dict[barcode]
            tag_length = len(barcode)+k
    return (region, tag_length)

class FailureTracker():
    '''
    Class to track failures at various stages in pipeline
    '''
    def __init__(self):
        self.fail_dict = {}

    def test(self, condition, description):

        # Add test to list of tests
        if not description in self.fail_dict:
            self.fail_dict[description] = 0

        # Increment fail_dict value
        if condition:
            self.fail_dict[description] += 1
            return True
        else:
            return False

    def get_results(self):
        return self.fail_dict
