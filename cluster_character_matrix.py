"""
Given a cell cluster it will give you the character matrix specific to that cluster
Input:
    1. Cell cluster Anndata as h5ad file
    2. All cell character matrix
Output:
    1. Modified Anndata ommiting the cells without any lineage data
    2. Character Matrix specific to that cluster
"""

# Environment settings
import scanpy as sc
import os
import pickle
import pandas as pd
import numpy as np


def ch_mat(adams, character_matrix):

    adams_modified = adams
    cluster_char_mat = []
    for i in range(adams.shape[0]):

        # for a cell find out the proper mutation state (one row of character matrix)
        cell = adams.obs_names.values[i]
        mutation_state = character_matrix[character_matrix.iloc[:, 1].isin([cell])]

        # sometimes we have two mutation states/cell, we just take the mutation associated with cell tumor
        if mutation_state.shape[0]>1:
            cell_tumor = adams.obs.at[cell, 'Tumor']
            mutation_state = mutation_state[mutation_state.iloc[:, 0] == cell_tumor]

            # if no tumor does not match then we took the cell as per tumor_ALL
            if mutation_state.shape[0] ==0:
                mutation_state = character_matrix[character_matrix.iloc[:, 1].isin([cell])][0:]
            cluster_char_mat.append(mutation_state)

        # if there is no lineage data, remove that cell
        elif mutation_state.shape[0]==0:
            adams_modified = adams_modified[adams_modified.obs.index != cell]

        else :
            cluster_char_mat.append(mutation_state)
    print()

    # we are making a final character matrix with all cells in the cluster
    cluster_char_mat_pd = pd.DataFrame()
    for ii in range(len(cluster_char_mat)):

        if ii == 0:
            row = cluster_char_mat[ii]

            # sometimes there are two mutation states with similar cell, we just took the first one
            # (AT) this may need special attention
            if row.shape[0]>1:
                cluster_char_mat_pd = row.iloc[[0]]
            else:
                cluster_char_mat_pd = row

        else:
            row = cluster_char_mat[ii]

            if row.shape[0]>1:
                cluster_char_mat_pd = (pd.concat
            ([cluster_char_mat_pd, row.iloc[[0]]], axis=0))
            else:
                cluster_char_mat_pd = (pd.concat
            ([cluster_char_mat_pd, row], axis=0))

    cluster_char_mat_pd = cluster_char_mat_pd.drop([0, 1], axis=1)
    cluster_char_mat_pd[cluster_char_mat_pd == -1] = np.nan
    return adams_modified, cluster_char_mat_pd

    print()



""""  ###### VARIABLES ########  """

main_dir = '../data/KPTracer-Data_divided/'
h5ad_data_path = main_dir + 'h5ad_files/yang_Early gastric_3457_Apc.h5ad'
character_matrix_path = main_dir + 'character_matrix/ALL_CELL_character_matrix_2.pkl'

""""  ###### VARIABLES ########  """

# read the data and character matrix
adams_old = sc.read(h5ad_data_path)
with open(character_matrix_path, 'rb') as file:
    character_matrix_old = pickle.load(file)


# process the character matrix/per cluster
adams, character_matrix \
    = ch_mat(adams_old, character_matrix_old)

main_save_dir = '../data/KPTracer-Data_divided/h5ad_filesWLineage/'
anndata_file_name = main_save_dir+'NEW_'+ h5ad_data_path.split('/')[-1]
CM_file_name = main_save_dir+'NEW_'+ (h5ad_data_path.split('/')[-1]).split('.')[0]+'_character_matrix.pkl'
adams.write_h5ad(anndata_file_name)
character_matrix.to_pickle(CM_file_name)
CSV_file_name = main_save_dir+'NEW_'+ (h5ad_data_path.split('/')[-1]).split('.')[0]+'_character_matrix.csv'
character_matrix.to_csv(CSV_file_name, index=False)
print()