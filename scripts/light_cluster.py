#!/usr/bin/env python3
"""
Corrects IGH only cloning with IGK/L annotations (VERSION 1)
"""

# Imports
import pandas as pd
import sys

# Parse arguments
heavy_file = sys.argv[1]
light_file = sys.argv[2]
cell_id = sys.argv[3]
clone_id = sys.argv[4]
out_file = sys.argv[5]


def clusterLinkage(cell_series, group_series):
    """
    Returns a dictionary of {cell_id : clone_id} that identifies clusters of cells by analyzing their 
    features (group_series).

    Arguments:
      cell_series : iter of cell_id's.
      group_series : iter of group_id's.

    Returns:
      assign_dict :  dictionary of {cell_id : clone_id}.
    """

    # assign initial clusters
    lc_dict = {}
    for cell, group in zip(cell_series, group_series):
        try:    
            lc_dict[group].append(cell)
        except KeyError:
            lc_dict[group] = [cell]

    # link clusters (ON^2) ...ie for cells with multiple light chains does complete linkage
    cluster_dict = {}
    for i, gene in enumerate(lc_dict.keys()):
        notfound = True
        for cluster in cluster_dict:
            if any(i in lc_dict[gene] for i in cluster_dict[cluster]):
                cluster_dict[cluster] = cluster_dict[cluster] + lc_dict[gene]
                notfound = False
                break
        if notfound:
            cluster_dict[i] = lc_dict[gene]
    
    assign_dict = {cell:k for k,v in cluster_dict.items() for cell in set(v)}
    
    return assign_dict


# read in heavy and light DFs
heavy_df = pd.read_csv(heavy_file, dtype = 'object', sep = '\t')
light_df = pd.read_csv(light_file, dtype = 'object', sep = '\t')


# generate a CELL:CLONE dictionary from heavy df and add to light df (basically an inner join)
clone_dict = {v[cell_id]:v[clone_id] for k,v in heavy_df[[clone_id, cell_id]].T.to_dict().items()}

light_df = light_df.loc[light_df[cell_id].apply(lambda x: x in clone_dict.keys()),]

light_df[clone_id] = light_df.apply(lambda row: clone_dict[row[cell_id]], axis = 1)

# generate a "cluster_dict" of CELL:CLONE dictionary from light df 
cluster_dict = clusterLinkage(light_df[cell_id], light_df.apply(lambda row: \
								  row['V_CALL'].split(',')[0].split('*')[0] +\
                                  ',' + row['J_CALL'].split(',')[0].split('*')[0] + ',' + \
                                  str(len(row['JUNCTION'])) + ',' + row[clone_id], axis = 1))


# add assignments to heavy_df
heavy_df = heavy_df.loc[heavy_df[cell_id].apply(lambda x: x in cluster_dict.keys()),:]

heavy_df[clone_id] = heavy_df[clone_id] + '_' + heavy_df.apply(lambda row: str(cluster_dict[row[cell_id]]), axis =1)

heavy_df.to_csv(out_file, sep = '\t')