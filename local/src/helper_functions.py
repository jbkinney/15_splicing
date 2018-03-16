import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
#import seaborn as sns
from scipy import stats

# Converts mm to inches
def mm2inch(value): 
    return value/25.4

# Converts width in mm to fraction of figure width
def width_mm2fig(width_mm,fig):
    fig_width_mm = fig.get_size_inches()[0]*25.4
    return width_mm/fig_width_mm

# Converts height in mm to fraction of figure height
def height_mm2fig(height_mm,fig):
    fig_height_mm = fig.get_size_inches()[1]*25.4
    return height_mm/fig_height_mm

# Class in which to store nodes for hierarchical annotation
class Node:
    def __init__(self, label, level, first, last, offset, children=[]):
        self.label = label
        self.level = level
        self.first = first
        self.last = last
        self.offset = offset
        self.children = children
        
    def __str__(self):
        s = '\t'*self.level + '%s [%s,%s]'%(self.label,self.first+self.offset,self.last+self.offset)
        return s
    
# Parse a dataframe into a hierarchical network of nodes
def parse_dataframe(df,level=0,offset=0):
    node_list = []
    column = df.iloc[:,0]
    for m in range(len(column)):
        entry = df.iloc[m,0]
        if m == 0:
            # Create new node and reset last_entry
            node_list.append(Node(label=entry,level=level,first=m,last=m,offset=offset))
            last_entry = entry
            
        elif entry == last_entry:
            # Adjust last attribute of most recent node
            node_list[-1].last = m
            
        else:
            # Compute children
            if df.shape[1]>1:
                first = node_list[-1].first
                last = node_list[-1].last
                new_df = df.iloc[first:(last+1),1:]
                node_list[-1].children = parse_dataframe(new_df,level=level+1,offset=offset+first)
                
            # Create new node and reset last_entry
            node_list.append(Node(label=entry,level=level,first=m,last=m,offset=offset))
            last_entry = entry
    
    # Get children of last node
    if df.shape[1]>1:
        first = node_list[-1].first
        last = node_list[-1].last
        new_df = df.iloc[first:(last+1),1:]
        node_list[-1].children = parse_dataframe(new_df,level=level+1,offset=offset+first)
            
    # Return node list
    return node_list

# Collect list of nodes
def get_nodes(top_node_list):
    all_nodes = []
    for node in top_node_list:
        all_nodes.append(node)
        all_nodes.extend(get_nodes(node.children))
    return all_nodes

# Ok, so here is a function to annotate the x-axis like a gel
def gelx(ax,df_annotation,annotation_spacing,fontsize=8,rotation=0):
    # Remove xticks
    ax.set_xticks([])
    
    # Compute annotation tree
    tree = parse_dataframe(df_annotation)

    # Get all nodes from tree
    all_nodes = get_nodes(tree)
    
    # Make labels
    num_levels = df_annotation.shape[1]
    num_rows = ax.get_ylim()[1]
    for node in all_nodes:
        y = -(num_rows + (num_levels - node.level - 1)*annotation_spacing)
        xa = node.first + node.offset + .2
        xb = node.last + node.offset + .8
        xc = 0.5*(xa+xb)
        plt.text(xc,y-.4*annotation_spacing,node.label,ha='center',va='center',fontsize=fontsize,rotation=rotation)
        if len(node.children) > 1 and len(node.label) > 0:
            plt.plot([xa,xb],[y,y],color='black',linewidth=1,clip_on=False)
            
    # Return axis
    return ax


# Ok, so here is a function to annotate the x-axis like a gel
def gely(ax,df_annotation,annotation_spacing,fontsize=8,\
    rotation=90,ha='center'):
    #if rotation == 0:
    #    ha='right'
    #else:
    #    ha='center'
        
    # Remove xticks
    ax.set_yticks([])
    
    # Compute annotation tree
    tree = parse_dataframe(df_annotation)

    # Get all nodes from tree
    all_nodes = get_nodes(tree)
    
    # Make labels
    num_levels = df_annotation.shape[1]
    num_rows = ax.get_ylim()[1]
    for node in all_nodes:
        x = (node.level - num_levels + 1)*annotation_spacing
        ya = (node.first + node.offset + .1) - num_rows
        yb = (node.last + node.offset + .9) - num_rows
        yc = 0.5*(ya+yb)
        plt.text(x-.4*annotation_spacing,yc,node.label,ha=ha,va='center',fontsize=fontsize,rotation=rotation)
        if len(node.children) > 1 and len(node.label) > 0:
            plt.plot([x,x],[ya,yb],color='black',linewidth=1,clip_on=False)
            
    # Return axis
    return ax