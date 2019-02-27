#!/usr/bin/env python3
"""
Corrects IGH only cloning with IGK/L annotations (VERSION 2)
"""

# Imports
import pandas as pd

# Parse arguments
heavy_file = sys.argv[1]
light_file = sys.argv[2]
cell_id = sys.argv[3]
clone_id = sys.argv[4]
out_file = sys.argv[5]


# read in heavy and light DFs
heavy_df = pd.read_csv(heavy_file, dtype = 'object', sep = '\t')
light_df = pd.read_csv(light_file, dtype = 'object', sep = '\t')

# add the VJJL column to light_df
light_df['VJJL'] = light_df.apply(lambda row: row['V_CALL'].split(',')[0].split('*')[0] +\
                                  ',' + row['J_CALL'].split(',')[0].split('*')[0] + ',' + \
                                  str(len(row['JUNCTION'])), axis = 1)

# generate a CELL:CLONE dictionary from heavy df and add to light df (basically an inner join)
clone_dict = {v[cell_id]:v[clone_id] for k,v in heavy_df[[clone_id, cell_id]].T.to_dict().items()}

light_df = light_df.loc[light_df[cell_id].apply(lambda x: x in clone_dict.keys()),]

light_df[clone_id] = light_df.apply(lambda row: clone_dict[row[cell_id]], axis = 1)

# identify the majority VJJL for each clone CLONE:VJJL
majority_cluster_dict = light_df.groupby(clone_id)\
    .apply(lambda row: max(list(row['VJJL']), key = list(row['VJJL']).count)).to_dict()

# assign each cell to the majority VJJL, else, use the standard VJJL to create CELL:VJJL
cluster_dict = light_df\
    .groupby(cell_id)\
    .apply(lambda row: majority_cluster_dict[list(row[clone_id])[0]] \
           if any(x == majority_cluster_dict[list(row[clone_id])[0]] for x in row['VJJL'])\
           else ";".join(row['VJJL'])).to_dict()

# add assignments to the heavy_df
heavy_df = heavy_df.loc[heavy_df[cell_id].apply(lambda x: x in cluster_dict.keys()),].reset_index()

heavy_df[clone_id] = heavy_df[clone_id] + '_' + heavy_df.apply(lambda row: cluster_dict[row[cell_id]], axis =1)

heavy_df.to_csv(out_file, sep = '\t')