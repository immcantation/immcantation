#!/usr/bin/env python3
"""
Corrects IGH only cloning with IGK/L annotations (VERSION 1)
"""

# Imports
import pandas as pd
import sys

# Parse arguments
#   1: heavy data
#   2: light data
#   3: output
#   4: TSV format. one of changeo (default if unspecified) or airr.
heavy_file = sys.argv[1]
light_file = sys.argv[2]
out_file = sys.argv[3]
try:
    format = sys.argv[4]
except IndexError:
    format = 'changeo'

# Set column names
if format == 'changeo':
    cell_id = 'CELL'
    clone_id = 'CLONE'
    v_call = 'V_CALL'
    j_call = 'J_CALL'
    junction = 'JUNCTION'
elif format == 'airr':
    cell_id = 'cell_id'
    clone_id = 'clone_id'
    v_call = 'v_call'
    j_call = 'j_call'
    junction = 'junction'
else:
    sys.exit("Invalid format %s" % format)


def clusterLinkage(cell_series, group_series):
    """
    Returns a dictionary of {cell_id : clone_id} that identifies clusters of cells by analyzing their shared
    features (group_series) using single linkage. 

    Arguments:
      cell_series : iter of cell_id's.
      group_series : iter of group_id's.

    Returns:
      assign_dict :  dictionary of {cell_id : clone_id}.
    """

    # assign initial clusters
    # initial_dict = {group1: [cell1, cell2]}
    initial_dict = {}
    for cell, group in zip(cell_series, group_series):
        try:    
            initial_dict[group].append(cell)
        except KeyError:
            initial_dict[group] = [cell]
               
    # single linkage clusters (ON^2) ...ie for cells with multiple light chains
    # cluster_dict = {1: [cell1, cell2]}, 2 cells belong in same group if they share 1 light chain 
    cluster_dict = {}
    for i, group in enumerate(initial_dict.keys()):
        cluster_dict[i] = initial_dict[group]
        for cluster in cluster_dict:
            # if initial_dict[group] and cluster_dict[cluster] share common cells, add initial_dict[group] to cluster
            if cluster != i and any(cell in initial_dict[group] for cell in cluster_dict[cluster]):
                cluster_dict[cluster] = cluster_dict[cluster] + initial_dict[group]
                del cluster_dict[i]
                break
    
    # invert cluster_dict for return
    assign_dict = {cell:k for k,v in cluster_dict.items() for cell in set(v)}
    
    return assign_dict


# read in heavy and light DFs
heavy_df = pd.read_csv(heavy_file, dtype = 'object', sep = '\t')
light_df = pd.read_csv(light_file, dtype = 'object', sep = '\t')

# filter multiple heavy chains
heavy_df = heavy_df.loc[heavy_df.groupby(cell_id)[cell_id].transform("count") == 1]

# transfer clone IDs from heavy chain df to light chain df
clone_dict = {v[cell_id]:v[clone_id] for k, v in heavy_df[[clone_id, cell_id]].T.to_dict().items()}
light_df = light_df.loc[light_df[cell_id].apply(lambda x: x in clone_dict.keys()), ]
light_df[clone_id] = light_df.apply(lambda row: clone_dict[row[cell_id]], axis = 1)

# generate a "cluster_dict" of CELL:CLONE dictionary from light df  (TODO: use receptor object V/J gene names)
cluster_dict = clusterLinkage(light_df[cell_id],
                              light_df.apply(lambda row: row[v_call].split(',')[0].split('*')[0] + \
                                             ',' + row[j_call].split(',')[0].split('*')[0] + ',' + \
                                             str(len(row[junction])) + ',' + row[clone_id], axis=1))

# add assignments to heavy_df
heavy_df = heavy_df.loc[heavy_df[cell_id].apply(lambda x: x in cluster_dict.keys()), :]
heavy_df[clone_id] = heavy_df[clone_id] + '_' + heavy_df.apply(lambda row: str(cluster_dict[row[cell_id]]), axis=1)

# write heavy chains
heavy_df.to_csv(out_file, sep='\t', index=False)