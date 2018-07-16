#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.image as mpimg

import os
import sys
import re
import pdb
import glob

print('Running plot_junctions.py ...')

sys.path.append('src')

# Specify input and output files
in_file = 'output/junction_displacements.txt'
out_files_dict = {
    'brca2':  'plots/junctions_brca2.pdf',
    'smn1':   'plots/junctions_smn1.pdf',
    'ikbkap': 'plots/junctions_ikbkap.pdf'
}

# For plotting
sns.set_style("whitegrid", {'axes.grid' : False})

# Load junction data
displacements_df = pd.read_csv(in_file,sep='\t')
displacements_df.set_index('experiment',drop=True,inplace=True)

# Get displacements values, as well as min and max
displacement_values = [int(x) for x in displacements_df.columns[1:-1]]
min_displacement = displacement_values[0]
max_displacement = displacement_values[-1]
xticks = [min_displacement-1] + \
         displacement_values + \
         [max_displacement+1]
xticklabels = ['<'] + list(displacements_df.columns[1:-1]) + ['>']
zero_index = xticks.index(0)

def plot_junction_data(subplot_loc_dict, subplot_dims, \
                       file_name, figsize):
    # Histogram junction lengths
    plt.figure(figsize=figsize)
    num_rows = subplot_dims[0]
    num_cols = subplot_dims[1]

    for experiment, loc in subplot_loc_dict.items():
        assert 0 < loc <= num_rows*num_cols,\
            'Error: subplot loc=%d invalid if num_rows=%d and num_cols=%d.'%\
            (loc,num_rows,num_cols)
        plt.subplot(num_rows,num_cols,loc)
        displacements = displacements_df.loc[experiment,:].values
        displacements = displacements/sum(displacements)
        plt.bar(xticks,displacements)
        plt.xticks(xticks,xticklabels)
        plt.yticks([0,.5,1])
        plt.ylim([0,1])
        plt.title('%s\n(%.0f%% correct)'%(experiment,100*displacements[zero_index]))
    plt.tight_layout()
    plt.savefig(file_name)
    print('Plot saved to %s.'%file_name)

# Plot results for brca2
subplot_loc_dict = {
'brca2_9nt_lib1_rep1':1,
'brca2_9nt_lib2_rep1':2,
}
plot_junction_data(subplot_loc_dict, 
    subplot_dims=(2,3), 
    figsize=(8,4), 
    file_name=out_files_dict['brca2'])

# Plot results for ikbkap
subplot_loc_dict = {
'ikbkap_9nt_lib1_rep1':1,
'ikbkap_9nt_lib1_rep2':2,
'ikbkap_9nt_lib1_rep3':3,
'ikbkap_9nt_lib2_rep1':4,
'ikbkap_9nt_lib2_rep2':5,
'ikbkap_9nt_lib2_rep3':6
}
plot_junction_data(subplot_loc_dict, 
    subplot_dims=(2,3), 
    figsize=(8,4), 
    file_name=out_files_dict['ikbkap'])


# Plot results for smn1
subplot_loc_dict = {
'smn1_9nt_lib1_rep1':1,
'smn1_9nt_lib1_rep2':2,
'smn1_9nt_lib1_rep3':3,
'smn1_9nt_lib2_rep1':4,
'smn1_9nt_lib2_rep2':5,
'smn1_9nt_lib2_rep3':6,
'smn1_9nt_lib3_rep1':7,
'smn1_9nt_lib3_rep2':8,
'smn1_9nt_lib3_rep3':9
}
plot_junction_data(subplot_loc_dict, 
    subplot_dims=(3,3), 
    figsize=(8,6), 
    file_name=out_files_dict['smn1'])

print('Done!')