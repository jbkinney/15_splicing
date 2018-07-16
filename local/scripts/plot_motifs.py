#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import os
import sys
import re
import pdb

print('Running plot_motifs.py ...')

sys.path.append('src')
import anylogo
import utils

# Converts splice site sequences to matrix form
bases = 'ACGU'
bases_mat = np.array(list(bases))


def seq_to_mat(seq):
    tmp_mat = np.array(list(seq))
    return (tmp_mat[:, np.newaxis] == bases_mat[np.newaxis, :]).astype(int)


# Load ratios
in_files_dict = {
    'lib_9nt': 'output/ratios_9nt_ss_lib.txt',
    'locus_9nt': 'output/ratios_9nt_ss_locus.txt',
}
in_human_file = 'data/human_splice_sites.xlsx'

matrix_dict = {}
psi_threshold = 20
cons_seq = 'CAGGUAAGU'

print('Loading human splice sites...')
# Load human splice sites
tmp_df = pd.read_excel(in_human_file, header=0)

# Extract GT splice sites
human_GU_df = tmp_df.iloc[:, 0:2]
human_GU_df.columns = ['human_ct', 'ss']
human_GU_df.dropna(how='any', axis=0, inplace=True)
human_GU_df['ss'] = [ss.replace('T', 'U') for ss in human_GU_df['ss']]
human_GU_df.set_index('ss', drop=True, inplace=True)

# Extract GC splice sites
human_GC_df = tmp_df.iloc[:, 6:8]
human_GC_df.columns = ['human_ct', 'ss']
human_GC_df.dropna(how='any', axis=0, inplace=True)
human_GC_df['ss'] = [ss.replace('T', 'U') for ss in human_GC_df['ss']]
human_GC_df.set_index('ss', drop=True, inplace=True)

# Concatenate GT and GC splice site data frames
human_df = pd.concat([human_GU_df, human_GC_df], axis=0, join='outer')
human_df['ss'] = [ss.replace('T', 'U') for ss in human_df.index]
human_df.set_index('ss', inplace=True, drop=True)

# Remove cryptic splice sites
human_df = utils.remove_cryptic_splice_sites(human_df)

# Compute matrices for human dfs
ss_len = len(human_df.index[0])
human_df_dict = {'all': human_df, 'GC': human_GC_df, 'GU': human_GU_df}
col_dict = {}
for name, df in human_df_dict.items():
    ct_mat = np.zeros([ss_len, len(bases)])
    for ss, row in df.iterrows():
        ct = row['human_ct']
        ss_mat = seq_to_mat(ss)
        ct_mat += ss_mat * ct
    col_dict[name] = ct_mat
matrix_dict['human'] = col_dict

# Comptue matrices for measured active sites
print('Running plot_motifs.py ...')
for name, in_file in in_files_dict.items():
    print('Loading %s...' % in_file)

    df = pd.read_csv(in_file, sep='\t')
    df.set_index('ss', inplace=True, drop=True)

    # Remove cryptic splice sites
    df = utils.remove_cryptic_splice_sites(df)

    site_length = len(df.index[0])

    # Process each column
    col_dict = {}
    for col in df.columns:
        tmp_df = df[[col]].dropna()

        # Normalize by consensus sequence
        if site_length == 9:
            assert cons_seq in tmp_df.index
            cons_seq_val = tmp_df.loc[cons_seq, col]

        else:
            raise Exception

        tmp_df.loc[:, col] = 100 * tmp_df.loc[:, col] / cons_seq_val

        # Get list of sites above threshold
        indices = tmp_df[col] > psi_threshold
        sites = tmp_df.index[indices].tolist()
        if len(sites) == 0:
            pdb.set_trace()
        assert len(sites) > 0
        assert len(sites[0]) == site_length

        # Fill counts matrix
        counts_matrix = np.zeros([site_length, 4])
        for ss in sites:
            counts_matrix += seq_to_mat(ss)

        col_dict[col] = counts_matrix
    matrix_dict[name] = col_dict

# Define plots to make:
plots_dict = {
    'plot1_all_9nt_loci': [
        ('locus_9nt', 'brca2_9nt', 'BRCA2'),
        ('locus_9nt', 'smn1_9nt', 'SMN1'),
        ('locus_9nt', 'ikbkap_9nt', 'IKBKAP'),
    ],
    'plot2_brca2_9nt_libs': [
        ('locus_9nt', 'brca2_9nt', 'BRCA2 all'),
        ('lib_9nt', 'brca2_9nt_lib1', 'BRCA2 lib1'),
        ('lib_9nt', 'brca2_9nt_lib2', 'BRCA2 lib2'),
    ],
    'plot3_smn1_9nt_libs': [
        ('locus_9nt', 'smn1_9nt', 'SMN1 all'),
        ('lib_9nt', 'smn1_9nt_lib1', 'SMN1 lib1'),
        ('lib_9nt', 'smn1_9nt_lib2', 'SMN1 lib2'),
        ('lib_9nt', 'smn1_9nt_lib3', 'SMN1 lib3'),
    ],
    'plot4_ikbkap_9nt_libs': [
        ('locus_9nt', 'ikbkap_9nt', 'IKBKAP all'),
        ('lib_9nt', 'ikbkap_9nt_lib1', 'IKBKAP lib1'),
        ('lib_9nt', 'ikbkap_9nt_lib2', 'IKBKAP lib2'),
    ],
    'plot5_human_9nt': [
        ('human', 'all', 'Human'),
        ('human', 'GU', 'Human, GU only'),
        ('human', 'GC', 'Human, GC only')
    ]
}

titles = list(plots_dict.keys())
titles.sort()


def sign(x):
    return '+' if x > 0 else '-'


for title in titles:
    plots = plots_dict[title]
    num_motifs = len(plots)
    fig = plt.figure(figsize=[4, 1.5 * num_motifs])

    for n, (name, col, header) in enumerate(plots):
        # extract matrix
        col_dict = matrix_dict[name]
        matrix = col_dict[col]
        num_sites = matrix.sum(axis=1)[0]

        motif = pd.DataFrame(data=matrix, columns=list('ACGU'))
        motif.index.name = 'pos'
        site_length = motif.shape[0]
        pos_values = motif.index
        tmp = [-3, -2, -1] + list(range(1, site_length - 2))
        pos_labels = ['%s %d' % (sign(x), abs(x)) for x in tmp]

        ax = fig.add_subplot(num_motifs, 1, n + 1)
        anylogo.draw(ax, prob_df=motif, logo_type='probability',
                     use_transparency=True, floor_line_width=1,
                     stack_order='big_on_top',
                     color_scheme='classic')

        if 'human' in name:
            clean_title = '%s (%d sites)' % (header, num_sites)
        else:
            clean_title = '%s, %d sites > %d PSI' % \
                      (header, num_sites, psi_threshold)
        plt.title(clean_title)
        plt.xticks(pos_values, pos_labels)
        plt.xlabel('')
        plt.yticks([])
        plt.ylabel('')

    plt.tight_layout()
    out_file = 'plots/motifs.%s.pdf' % title
    plt.savefig(out_file)
    print('Output saved to %s' % out_file)

print('Done!')
