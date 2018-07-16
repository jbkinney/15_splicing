#!/usr/bin/env python
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import os
import sys
import re
import pdb
import glob
from scipy.special import logit

print('Running plot_correlations.py ...')

sys.path.append('src')
import utils
from helper_functions import gelx, gely

os.system('rm plots/correlations*.pdf')

def get_annotation_df(raw_names, key):
    if ('rep' in key) or ('all' in key):
        annotation_df = pd.DataFrame(columns=['locus', 'lib', 'rep'])
        pattern = re.compile('(?P<locus>\S+)_.*nt' +
                             '_lib(?P<lib_num>[0-9])' +
                             '_rep(?P<rep_num>[0-9])')
        for n, name in enumerate(raw_names):
            m = pattern.match(name)
            if not m:
                raise Exception
                # clean_names.append(name)
                # continue
            locus = m.group('locus')
            if '11nt' in locus:
                locus += ' 11nt'
            lib = 'lib #%s' % str(m.group('lib_num'))
            rep = '%s' % str(m.group('rep_num'))
            locus = locus.upper() #locus[0].upper() + locus[1:]
            annotation_df.loc[n] = {'locus': locus, 'lib': lib, 'rep': rep}
            annotation_spacing = 1

    elif ('lib' in key):
        annotation_df = pd.DataFrame(columns=['locus', 'lib'])
        pattern = re.compile('(?P<locus>\S+)_.*nt' +
                             '_lib(?P<lib_num>[0-9])')
        for n, name in enumerate(raw_names):
            m = pattern.match(name)
            if not m:
                raise Exception
                # clean_names.append(name)
                # continue
            locus = m.group('locus')
            if '11nt' in locus:
                locus += ' 11nt'
            lib = 'lib #%s' % str(m.group('lib_num'))
            locus = locus.upper() #locus[0].upper() + locus[1:]
            annotation_df.loc[n] = {'locus': locus, 'lib': lib}
            annotation_spacing = .4

    elif ('locus' in key):
        annotation_df = pd.DataFrame(columns=['locus'])
        pattern = re.compile('(?P<locus>\S+)_.*nt')
        for n, name in enumerate(raw_names):
            m = pattern.match(name)
            if not m:
                locus = name
            else:
                locus = m.group('locus')
                if '11nt' in locus:
                    locus += ' 11nt'
                locus = locus.upper() #locus[0].upper() + locus[1:]
            annotation_df.loc[n] = {'locus': locus}
            annotation_spacing = .2
    else:
        assert False, 'Error! unparsable raw_names'

    return annotation_df, annotation_spacing


def clean_sample_names(raw_names, key):
    clean_names = []
    if ('rep' in key) or ('all' in key):
        pattern = re.compile('(?P<name>\S+)_.*nt' +
                             '_lib(?P<lib_num>[0-9])' +
                             '_rep(?P<rep_num>[0-9])')
        for name in raw_names:
            m = pattern.match(name)
            if not m:
                clean_names.append(name)
                continue
            name = m.group('name')
            if '11nt' in name:
                name += ' 11nt'
            lib_num = int(m.group('lib_num'))
            rep_num = int(m.group('rep_num'))
            #tmp = name[0].upper() + name[1:]
            tmp = name.upper()
            clean_names.append( \
                '%s (%d.%d)' % (tmp, lib_num, rep_num))
    elif ('lib' in key):
        pattern = re.compile('(?P<name>\S+)_.*nt' +
                             '_lib(?P<lib_num>[0-9])')
        for name in raw_names:
            m = pattern.match(name)
            if not m:
                clean_names.append(name)
                continue
            name = m.group('name')
            if '11nt' in name:
                name += ' 11nt'
            lib_num = int(m.group('lib_num'))
            #tmp = name[0].upper() + name[1:]
            tmp = name.upper()
            clean_names.append( \
                '%s (%d)' % (tmp, lib_num))
    elif ('locus' in key):
        pattern = re.compile('(?P<name>\S+)_.*nt')
        for name in raw_names:
            m = pattern.match(name)
            if not m:
                clean_names.append(name)
                continue
            name = m.group('name')
            if '11nt' in name:
                name += ' 11nt'
            #tmp = name[0].upper() + name[1:]
            tmp = name.upper()
            clean_names.append('%s' % (tmp))
    return clean_names


in_files_dict = {
    'all_9nt': 'output/ratios_9nt_ss_all.txt',
    'rep_9nt': 'output/ratios_9nt_ss_rep.txt',
    'lib_9nt': 'output/ratios_9nt_ss_lib.txt',
    'locus_9nt': 'output/ratios_9nt_ss_locus.txt',
}
description_dict = {
    'all_9nt': 'all 9nt replicates',
    'rep_9nt': 'selected 9nt replicates',
    'lib_9nt': '9nt libraries',
    'locus_9nt': '9nt loci',
}

# Set style
# sns.set_style("whitegrid", {'axes.grid' : False})
fontsize = 9
cmap = cm.get_cmap('Spectral')
# cmap = sns.cubehelix_palette(8, start=1.2, rot=0.0, reverse=True, as_cmap=True)

# Consider 3 different subsets of sites
for subset in [None]:  # , 'GU', 'GC']:

    # Compute correlations for each data set
    for in_key, in_file_name in in_files_dict.items():
        ratios_df = pd.read_csv(in_file_name, sep='\t')
        ratios_df.set_index('ss', inplace=True, drop=True)

        # If restricting to a subset of sites
        if subset:
            indices = [x[3:5] == subset for x in ratios_df.index]
            ratios_df = ratios_df.loc[indices, :]

        # Remove cryptic sites
        ratios_df = utils.remove_cryptic_splice_sites(ratios_df)

        # Make title
        title = '$R^2$ between %s' % description_dict[in_key]
        if subset:
            title += ' %s only' % subset

        # Set out_file name
        if subset:
            out_file_name = 'plots/correlations_%s_%s.pdf' % (in_key, subset)
        else:
            out_file_name = 'plots/correlations_%s.pdf' % (in_key)

        # Reorder columns
        cols = ratios_df.columns
        new_cols = \
            [c for c in cols if 'brca2' in c] + \
            [c for c in cols if 'smn1' in c] + \
            [c for c in cols if 'ikbkap' in c]
        ratios_df = ratios_df[new_cols]

        # Clean sample names
        sample_names = \
            clean_sample_names(ratios_df.columns, in_key)
        num_samples = len(sample_names)
        annotation_df, annotation_spacing = \
            get_annotation_df(ratios_df.columns, in_key)
        if num_samples == 3 and annotation_df.shape[1] > 1:
            annotation_spacing = .2  # This is a complete hack
        if num_samples == 9 and annotation_df.shape[1] > 1:
            annotation_spacing = .5  # This is a complete hack

        # Comptue correlations
        corr_mat = ratios_df.corr(method='pearson') ** 2

        # Draw heatmap
        fig = plt.figure(figsize=(6, 5))
        left = .15
        bottom = .05
        width = .8
        height = .8
        ax = fig.add_axes([left, bottom, width, height],aspect='equal')
        sns.heatmap(
            corr_mat, annot_kws={"size": 7},
            cmap=cmap, cbar_kws={"pad": .05, "aspect": 50},
            vmin=0, vmax=1)
        gelx(ax, annotation_df,
             fontsize=fontsize,
             annotation_spacing=annotation_spacing)
        gely(ax, annotation_df,
             annotation_spacing=annotation_spacing,
             fontsize=fontsize)

        # Change fontsize of colorbar tick labels
        cax = plt.gcf().axes[-1]
        cax.tick_params(labelsize=fontsize)

        # Display correlation numbers
        if num_samples < 8:
            for i in range(num_samples):
                for j in range(num_samples):
                    corr = corr_mat.iloc[i, j]
                    if .4 < corr < .6:
                        color = 'dimgray'
                    else:
                        color = 'white'
                    # plt.text(i + .5, num_samples - j - .5, '%0.2f' % corr,
                    #          color=color, fontsize=fontsize,
                    #          horizontalalignment='center',
                    #          verticalalignment='center')
                    plt.text(i + .5, j + .5, '%0.2f' % corr,
                             color=color, fontsize=fontsize,
                             horizontalalignment='center',
                             verticalalignment='center')

        # Save outfile
        plt.savefig(out_file_name)

        # plt.close()
        ax.set_title(title, fontsize=fontsize)
        ttl = ax.title
        ttl_y_pos = 1.03+.04*annotation_df.shape[1]
        ttl.set_position([.5, ttl_y_pos])
        print('Plot written to %s' % out_file_name)

# Load computational predictions
comp_file = 'data/scored_splice_sites.txt'
comp_df = pd.read_csv(comp_file,
                      usecols=[0, 2, 4, 6, 8],
                      names=['ss', 'MaxEnt', 'MDD', 'MM', 'WMM'],
                      delim_whitespace=True)
comp_df['ss'] = [ss.replace('T', 'U') for ss in comp_df['ss']]
comp_df.set_index('ss', drop=True, inplace=True)

# Add in RNAhybrid predictions
rnahybrid_file = 'output/rnahybrid_predictions.txt'
tmp_df = pd.read_csv(rnahybrid_file, sep='\t', names=['ss', 'RNAhyb'])
tmp_df.set_index('ss', drop=True, inplace=True)
comp_df = comp_df.merge(tmp_df, left_index=True, right_index=True, how='outer')

# Smush predictions
comp_df = utils.smush_predictions(comp_df)

# Remove cryptic sites
comp_df = utils.remove_cryptic_splice_sites(comp_df)

comp_methods = ['all', 'MaxEnt', 'MDD', 'MM', 'WMM', 'RNAhyb']
for comp_method in comp_methods:

    # If only doing comparison to maxent
    if not comp_method == 'all':
        tmp_df = comp_df[[comp_method]].copy()
    else:
        tmp_df = comp_df.copy()

    # Load locus correlations, with only GT sites
    in_key = 'locus_9nt'
    in_file_name = in_files_dict[in_key]
    ratios_df = pd.read_csv(in_file_name, sep='\t')
    ratios_df.set_index('ss', drop=True, inplace=True)

    # Remove cryptic splice sites
    ratios_df = utils.remove_cryptic_splice_sites(ratios_df)

    # Reorder columns
    cols = ratios_df.columns
    new_cols = \
        [c for c in cols if 'brca2' in c] + \
        [c for c in cols if 'smn1' in c] + \
        [c for c in cols if 'ikbkap' in c]
    ratios_df = ratios_df[new_cols]

    # Merge
    both_df = ratios_df.merge(tmp_df, left_index=True, right_index=True, how='left')

    # Restrict to GT splice sites
    indices = [x[3:5] == 'GU' for x in both_df.index]
    both_df = both_df.loc[indices, :]

    # Clean sample names
    sample_names = \
        clean_sample_names(both_df.columns, in_key)
    num_samples = len(sample_names)
    annotation_df, annotation_spacing = \
        get_annotation_df(both_df.columns, in_key)
    if len(both_df.columns) > 5:
        annotation_spacing = .5  # This is a complete hack

    # Make title
    title = r'$R^2$, measurements & predictions'

    # Set out_file name
    out_file_name = 'plots/correlations_%s_comp_%s.pdf' % (in_key, comp_method.replace(' ', '_'))

    # Comptue correlations
    corr_mat = both_df.corr(method='pearson') ** 2

    # Draw heatmap
    fig = plt.figure(figsize=(6, 5))
    left = .15
    bottom = .05
    width = .8
    height = .8
    ax = fig.add_axes([left, bottom, width, height], aspect='equal')
    sns.heatmap(
        corr_mat, annot_kws={"size": 7},
        cmap=cmap, cbar_kws={"pad": .05, "aspect": 50},
        vmin=0, vmax=1)
    gelx(ax, annotation_df,
         annotation_spacing=annotation_spacing,
         fontsize=fontsize)
    gely(ax, annotation_df,
         annotation_spacing=annotation_spacing,
         fontsize=fontsize)

    # Change fontsize of colorbar tick labels
    cax = plt.gcf().axes[-1]
    cax.tick_params(labelsize=fontsize)

    # Display correlation numbers
    if num_samples < 8:
        for i in range(num_samples):
            for j in range(num_samples):
                corr = corr_mat.iloc[i, j]
                if .4 < corr < .6:
                    color = 'dimgray'
                else:
                    color = 'white'
                # plt.text(i + .5, num_samples - j - .5, '%0.2f' % corr,
                #          color=color, fontsize=fontsize,
                #          horizontalalignment='center',
                #          verticalalignment='center')
                plt.text(i + .5, j + .5, '%0.2f' % corr,
                          color=color, fontsize=fontsize,
                          horizontalalignment='center',
                          verticalalignment='center')

    # Render title
    ax.set_title(title, fontsize=fontsize)
    ttl = ax.title
    ttl_y_pos = 1.03 + .04 * annotation_df.shape[1]
    ttl.set_position([.5, ttl_y_pos])

    # Save outfile
    plt.savefig(out_file_name)

    # plt.close()
    print('Plot written to %s' % out_file_name)

print('Done!')
