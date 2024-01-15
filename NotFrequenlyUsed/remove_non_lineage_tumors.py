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
character_matrix_path = '../../data/KPTracer-Data/trees/'

""""  ###### VARIABLES ########  """

tumor_no_lineage = []
file_path = "report2_refined.txt"

with open(file_path, 'r') as text_file:
    # Read the entire contents of the file into a string
    for line in text_file:
        tumor_no_lineage.append(line.replace(' \n', ''))

# for all interested cell clusters, find the unique tumor classification
all_files = os.listdir(h5ad_data_path)
files = [file for file in all_files if not file.startswith('.DS_Store')]

file_path = "report2_new.txt"
with open(file_path, 'w') as report_file:
    for file in files:

        unique_tumors = []
        print("cell cluster: ", file)
        report_file.write(f"\ncell cluster {file}")

        final_character_matrix = np.zeros((2, 2))
        lineage_not_found = []
        adams = sc.read(os.path.join(h5ad_data_path, file))
        unique_tumor_values = adams.obs['Tumor'].unique()
        tumor_num = 0

        for tumor in unique_tumor_values:
            unique_tumors.append(tumor)

        set1 = set(unique_tumors)
        set2 = set(tumor_no_lineage)
        common_items = set1.intersection(set2)
        common_list = list(common_items)

        if len(common_list):
            adams_new = adams
            for tumors in common_list:
                adams_new = adams_new[~adams_new.obs['Tumor'].isin([tumors])].copy()

            file = "../new_h5ad_data/new_" + file
            adams_new.write_h5ad(os.path.join(h5ad_data_path, file))
            report_file.write(f"\nold file cell num {adams.shape}")
            report_file.write(f"\nnew file cell num {adams_new.shape}")
            print(f"old file cell num {adams.shape}")
            print(f"new file cell num {adams_new.shape}")

    #
    #
    #
    #
    #
    #     print("num of unique tumor ", len(unique_tumor_values))
    #     report_file.write(f"\nnum of unique tumor {unique_tumor_values}")
    #
    #     # Yang provided the character matrix as per tumors.
    #     # for all tumor values, read the character matrix.
    #     for value in unique_tumor_values:
    #         tumor_character_matrix = []
    #         tumor_file_name = os.path.join(character_matrix_path, (value + "_character_matrix.txt"))
    #
    #         with open(os.path.join(tumor_file_name), 'r') as text_file:
    #             # Read the entire contents of the file into a string
    #             for line in text_file:
    #                 print()
    #
    #             with open(os.path.join(tumor_file_name), 'r') as text_file:
    #                 line_num = 0
    #
    #                 # each line for each cell
    #                 for line in text_file:
    #                     if line_num == 0:
    #                         print()
    #                     else:
    #                         line_component = line.strip().split('\t')
    #
    #                         # except the cell name make mutation score integer from str
    #                         for i in range(1, len(line_component)):
    #                             try:
    #                                 line_component[i] = float(line_component[i])
    #                             except:
    #                                 line_component[i] = -1
    #                         tumor_character_matrix.append(line_component)
    #                     line_num += 1
    #
    #             # final character matrix for all cell in that cluster
    #             tumor_character_matrix = np.array(tumor_character_matrix)
    #
    #             final_cell_clm_num = np.shape(final_character_matrix)[1]
    #             current_cell_clm_num = np.shape(tumor_character_matrix)[1]
    #
    #             if tumor_num > 0:
    #                 if (final_cell_clm_num > current_cell_clm_num):  # if previous tumor cells have more mutation state
    #                     matrix_padding = np.full((np.shape(tumor_character_matrix)[0],
    #                                               abs(final_cell_clm_num - current_cell_clm_num)), 0)
    #                     tumor_character_matrix = np.hstack((tumor_character_matrix, matrix_padding))
    #                 elif (final_cell_clm_num < current_cell_clm_num):  # if current tumor cells have more mutation state
    #                     matrix_padding = np.full((np.shape(final_cell_clm_num)[0],
    #                                               abs(final_cell_clm_num - current_cell_clm_num)), 0)
    #                     tumor_character_matrix = np.hstack((tumor_character_matrix, matrix_padding))
    #
    #                 final_character_matrix = np.concatenate((final_character_matrix, tumor_character_matrix), axis=0)
    #
    #             else:  # for first tumor appending
    #                 final_character_matrix = tumor_character_matrix
    #             print()
    #
    #             lineage_not_found.append(value)
    #         tumor_num += 1
    #
    #     final_character_matrix = pd.DataFrame(final_character_matrix)
    #     final_name = os.path.join(h5ad_data_path,
    #                               ("../character_matrix/character_matrix_" + file.split(".").pop(0)) + ".pkl")
    #     final_character_matrix.to_pickle(final_name)
    #
    #     if len(lineage_not_found):
    #         adams_new = adams
    #         for tumors in lineage_not_found:
    #             print("lineage data NOT FOUND ", tumors)
    #             report_file.write(f"\nlineage data NOT FOUND {tumors}")
    #
    #             adams_new = adams_new[~adams_new.obs['Tumor'].isin([tumors])].copy()
    #
    #         file = "../new_h5ad_data/new_" + file
    #         adams_new.write_h5ad(os.path.join(h5ad_data_path, file))
    #         report_file.write(f"\nold file cell num {adams.shape}")
    #         report_file.write(f"\nnew file cell num {adams_new.shape}")
    #         print(f"old file cell num {adams.shape}")
    #         print(f"new file cell num {adams_new.shape}")
    #
    #         print()
