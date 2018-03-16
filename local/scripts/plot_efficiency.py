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

print('Running plot_efficiency.py ...')

sys.path.append('src')

# Specify input and output files
in_file_dict = {
    'num_reads':'output/efficiency_by_stage_num_reads.txt',
    'filtering':'output/efficiency_by_stage_filtering.txt'
}
out_files_dict = {
    'num_reads':'plots/efficiency_num_reads.pdf',
    'filtering':'plots/efficiency_filtering.pdf',
}

# For plotting
sns.set_style("whitegrid", {'axes.grid' : False})

#
# Plot total number of reads
#

# Load efficiency data by stage
df = pd.read_csv(in_file_dict['num_reads'],sep='\t')

# Make figure
plt.figure(figsize=[8,5])
    
# Make barplot
ax = sns.barplot(x='label',y='num_reads', hue='stage', data=df)

# Fix x labels
xticklabels = list(df['label'][df['stage']=='raw'])
ax.set_xticklabels(xticklabels,rotation=-45,horizontalalignment='left')

# Add legend
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12), ncol=3)
ax.set_xlabel('Analysis group + LID')
ax.set_ylabel('Number of reads')
ax.set_title('Number of reads passing filters',y=1.1)
plt.tight_layout()
plt.savefig(out_files_dict['num_reads'])
print('Plot saved to %s.'%out_files_dict['num_reads'])


#
# Plot percentage of reads passing filters.
#

# Load efficiency data by stage
df = pd.read_csv(in_file_dict['filtering'],sep='\t')

# Make figure
plt.figure(figsize=[8,5])
    
# Make barplot
ax = sns.barplot(x='label',y='frac', hue='stage', data=df)

# Fix x labels
xticklabels = list(df['label'][df['stage']=='sorted'])
ax.set_xticklabels(xticklabels,rotation=-45,horizontalalignment='left')

# Add legend
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12), ncol=2)
ax.set_xlabel('Analysis group + LID')
ax.set_ylabel('Efficiency')
ax.set_title('Efficiency of filtering steps',y=1.1)
plt.tight_layout()
plt.savefig(out_files_dict['filtering'])
print('Plot saved to %s.'%out_files_dict['filtering'])

print('Done!')
