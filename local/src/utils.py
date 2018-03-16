from __future__ import division
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

import os
import sys
import re
import pdb
import glob

WORKING_DIR = '/Users/jkinney/Dropbox/15_splicing/pipeline_local/'

RESULT_FILE_PATTERN = \
    re.compile('(?P<directory>\S+)/results\.(?P<experiment>.*)\.txt')
EXPERIMENT_PATTERN = \
    re.compile('(?P<locus>.+)_lib(?P<lib_num>[0-9]+)_rep(?P<rep_num>[0-9]+)')
COUNT_FILE_PATTERN = \
    re.compile('(?P<directory>\S+)/counts\.(?P<sample>.*)\.txt')


class Result:
    ''' Class to hold information on a results file'''

    def __init__(self, file_name):
        self.valid = False
        m = re.match(RESULT_FILE_PATTERN, file_name)
        if m:
            self.data_file = file_name

            d1 = m.groupdict()
            self.directory = d1['directory']
            self.experiment = d1['experiment']
            self.report_file = '%s/report.results.%s.txt' % (self.directory, self.experiment)

            d2 = re.match(EXPERIMENT_PATTERN, self.experiment).groupdict()
            self.locus = d2['locus']
            self.lib_num = int(d2['lib_num'])
            self.rep_num = int(d2['rep_num'])
            self.valid = True

            # Read in data
            assert os.path.isfile(self.data_file)
            self.data_df = pd.read_csv(self.data_file, delim_whitespace=True)

            # Load report file
            assert os.path.isfile(self.report_file)
            with open(self.report_file) as f:
                self.report = f.read()

    def __repr__(self):
        return '%s_lib%d_rep%d' % (self.locus, self.lib_num, self.rep_num)


def fill_results_dict(files_glob):
    '''
    Loads results from BNB into a results dict 
    '''
    results_dict = {}
    file_names = glob.glob(files_glob)
    assert len(file_names) > 0, 'Error! No files found for glob %s.' % files_glob

    # Load files
    for file_name in file_names:
        result = Result(file_name)
        if result.valid:
            name = repr(result)
            results_dict[name] = result
            print('Loading %s ...' % name)
    print('Done!')

    return results_dict


def get_ss_counts(results_glob, ss_length):
    # Create sequences
    if ss_length == 9:
        ss_list = make_iupac_seqs('NNNGYNNNN')
    elif ss_length == 11:
        ss_list = make_iupac_seqs('NNNGYNNNNNN')
    else:
        print('Error! ss_length=%d is invalid' % ss_length)
        raise Exception

    # Fill results dict
    file_names = glob.glob(results_glob)
    assert len(file_names) > 0, \
        'Error! No files found for glob %s.' % files_glob

    # Load files
    counts_df = pd.DataFrame(index=ss_list)
    for file_name in file_names:
        # Create results object
        result = Result(file_name)
        assert result.valid, 'Error! Invalid result!'
        name = repr(result)
        print('get_ss_counts: %s loaded.' % name)

        # Marginalize over barcodes; keep only ex and tot cts. 
        x = result.data_df.groupby(['ss']).sum()

        # Fill in counts values
        counts_df[name + '_ex_ct'] = x['ex_ct'].astype(int)
        counts_df[name + '_tot_ct'] = x['tot_ct'].astype(int)
        counts_df[name + '_lib_ct'] = x['lib_ct'].astype(int)
        counts_df[name + '_mis_ct'] = x['mis_ct'].astype(int)

    # Fill NaN entries with 0 
    counts_df.fillna(0, inplace=True)

    # Remove unexpected splice sites
    counts_df = counts_df.loc[ss_list, :]

    # Transform all 'T' to 'U' in index
    counts_df['ss'] = [s.replace('T', 'U') for s in counts_df.index]
    counts_df.set_index('ss', drop=True, inplace=True)
    counts_df.sort_index()

    return counts_df


def splice(names, elements, unique=False):
    '''
    Splices either a single name or a list of names.
    Elements is a list of components to keep.
    '''
    if type(names) == str:
        output_list = False
        names = [names]
    else:
        output_list = True

    new_names = \
        ['_'.join([name.split('_')[x] for x in elements]) for name in names]
    if unique:
        new_names = list(set(new_names))
        new_names.sort()

    if output_list:
        return new_names
    else:
        return new_names[0]


def counts_to_ratios(counts_df, min_tot_ct=5):
    '''
    Given a counts data frame, returns a data frame of normalized
    experiment ratios
    '''

    # Get experiments
    cols = counts_df.columns
    experiments = set([])
    for col in cols:
        m = re.search(EXPERIMENT_PATTERN, col)
        experiments.add(m.group(0))
    experiments = list(experiments)
    experiments.sort()

    # Compute ratios of normalized counts
    ratios_df = pd.DataFrame(index=counts_df.index)
    for name in experiments:
        # Get counts for ex and tot samples
        ex_cts = counts_df[name + '_ex_ct'].fillna(0)
        tot_cts = counts_df[name + '_tot_ct'].fillna(0)

        # Compute total number of identified reads for ex and tot samples
        ex_ct_sum = sum(ex_cts)
        tot_ct_sum = sum(tot_cts)

        # Compute indices with valid min tot_ct valies
        indices = tot_cts >= min_tot_ct

        # Compute normalized ratios
        ratios_df[name] = ((ex_cts / tot_cts) * (tot_ct_sum / ex_ct_sum))[indices]

    # Sort columns
    cols = list(ratios_df.columns)
    cols.sort()
    ratios_df = ratios_df[cols]

    return ratios_df


def make_iupac_seqs(iupac_motif):
    """Returns a list of all DNA sequences that match an IUPAC motif"""
    iupac_dict = {
        'A': 'A',
        'C': 'C',
        'G': 'G',
        'T': 'T',
        'R': 'AG',
        'Y': 'CT',
        'S': 'GC',
        'W': 'AT',
        'K': 'GT',
        'M': 'AC',
        'B': 'CGT',
        'D': 'AGT',
        'H': 'ACT',
        'V': 'ACG',
        'N': 'ACGT'
    }

    if len(iupac_motif) == 0:
        return ['']
    else:
        # Recursively get subseqs
        suffix_seqs = make_iupac_seqs(iupac_motif[1:])

        # Get prefix characters
        key = iupac_motif[0]
        assert key in iupac_dict.keys(), 'Invalid character: %s' % key
        prefix_chars = iupac_dict[key]

        # Build new sequences
        seqs = [c + s for c in prefix_chars for s in suffix_seqs]
        seqs.sort()
        return seqs


def kraskov_MI(xs, ys, k_neighbors=6, noise=1E-8):
    N = len(xs)
    assert len(ys) == N
    tmp_df = pd.DataFrame()
    tmp_df['x'] = xs
    tmp_df['y'] = ys
    in_file_name = '.kraskov_input_df.txt.tmp'
    out_file_name = '.mi.txt.tmp'
    tmp_df.to_csv(in_file_name, sep='\t', header=False, index=False)

    command = """
    MIhigherdim %s 2 1 1 %d %d %f > %s
    """ % (in_file_name, N, k_neighbors, noise, out_file_name)
    os.system(command)

    f = open(out_file_name, 'r')
    mi = float(f.readline().strip())
    return mi


def smush_predictions(prediction_df):
    new_df = prediction_df.copy()
    cols = prediction_df.columns
    for col in cols:
        if col == 'MaxEnt':
            floor = 0
        elif col == 'MDD':
            floor = 5
        elif col == 'MM':
            floor = 3
        elif col == 'WMM':
            floor = 1
        else:
            floor = -np.Inf
        if floor > -np.Inf:
            new_df.loc[new_df[col] < floor, col] = floor
    return new_df


def remove_cryptic_splice_sites(in_df, column=None):
    if column is None:
        vals = in_df.index
    else:
        assert column in in_df.columns, 'Error! column %s not found in in_df.columns.' % column
        vals = in_df.columns

    indices = [(ss[1:3] != 'GU') and (ss[1:3] != 'GT') for ss in vals]
    out_df = in_df.loc[indices, :]
    return out_df
