"""
Given the .h5ad files it goes to Yang character matrices and form a separate character matrix file  for that particular cell clluster
"""

# Environment settings
import scanpy as sc
import os
import numpy as np
import pandas as pd

""""  ###### VARIABLES ########  """

h5ad_data_path = '../data/KPTracer-Data_divided/h5ad_files/'
character_matrix_path = '../data/KPTracer-Data/trees/'

""""  ###### VARIABLES ########  """

# for all interested cell clusters, find the unique tumor classification
all_files = os.listdir(h5ad_data_path)
files = [file for file in all_files if not file.startswith('.DS_Store')]


for file in files:
    #print("cell cluster: ", file)

    final_character_matrix = np.zeros((2,2))
    lineage_not_found = []
    adams = sc.read(os.path.join(h5ad_data_path, file))
    unique_tumor_values = adams.obs['Tumor'].unique()
    tumor_num =0
    #print("num of unique tumor ", len(unique_tumor_values))
    # Yang provided the character matrix as per tumors.
    # for all tumor values, read the character matrix.
    for value in unique_tumor_values:
        tumor_character_matrix = []
        tumor_character_matrix_cellspecific = []

        tumor_file_name = os.path.join(character_matrix_path, (value+"_character_matrix.txt"))
        exists = os.path.exists(tumor_file_name)
        if exists:
            print("put some thing")
        else:
            tumor_file_name = os.path.join(character_matrix_path,
                    (line.split("_").pop(0) + "_" + line.split("_").
                        pop(1) + "_All_character_matrix.txt"))

        # cell transcriptome as per thr tumor
        adam_cells = adams[adams.obs['Tumor'].isin([value])].obs_names

        with open(tumor_file_name, 'r') as text_file:
            line_num = 0
            # each line for each cell
            for line in text_file:
                if line_num == 0:
                    print()
                else:
                    line_component = line.strip().split('\t')
                    # except the cell name make mutation score integer from str
                    for i in range(1, len(line_component)):
                        try:
                            line_component[i] = float(line_component[i])
                        except:
                            line_component[i] = -1
                    tumor_character_matrix.append(line_component)
                line_num += 1

        # final character matrix for all cell in that cluster
        tumor_character_matrix = np.array(tumor_character_matrix)

        # just keep the cells we want from the tumor chracter matrix
        for cell in adam_cells:
            for tumor_cell in tumor_character_matrix:
                if cell == tumor_cell[0]:
                    tumor_character_matrix_cellspecific.append(tumor_cell)

        final_cell_clm_num = np.shape(final_character_matrix)[1]
        current_cell_clm_num = np.shape(tumor_character_matrix)[1]

        if tumor_num > 0:
            if (
                    final_cell_clm_num > current_cell_clm_num):  # if previous tumor cells have more mutation state
                matrix_padding = np.full((np.shape(tumor_character_matrix)[0],
                                          abs(final_cell_clm_num - current_cell_clm_num)), 0)
                tumor_character_matrix = np.hstack((tumor_character_matrix, matrix_padding))
            elif (
                    final_cell_clm_num < current_cell_clm_num):  # if current tumor cells have more mutation state
                matrix_padding = np.full((np.shape(final_cell_clm_num)[0],
                                          abs(final_cell_clm_num - current_cell_clm_num)), 0)
                tumor_character_matrix = np.hstack((tumor_character_matrix, matrix_padding))

            final_character_matrix = np.concatenate((final_character_matrix, tumor_character_matrix),
                                                    axis=0)
        else:  # for first tumor appending
            final_character_matrix = tumor_character_matrix

        tumor_num += 1

        final_character_matrix = pd.DataFrame(final_character_matrix)
        final_name = os.path.join(h5ad_data_path,
                                  ("../character_matrix/character_matrix_" + file.split(".").pop(
                                      0)) + ".pkl")
        final_character_matrix.to_pickle(final_name)
        print()

