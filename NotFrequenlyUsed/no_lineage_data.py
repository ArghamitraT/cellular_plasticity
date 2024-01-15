import os

character_matrix_path = '../../data/KPTracer-Data/trees/'

file_path = "report2_refined.txt"
with open(file_path, 'w') as report_file:
    with open('report2.txt', 'r') as text_file:
        # Read the entire contents of the file into a string
        for line in text_file:
            try:  # (some tumor cells do not have lineage data)
                with open(os.path.join(character_matrix_path,
                (line.split("_").pop(0) + "_" + line.split("_").
                        pop(1) + "_All_character_matrix.txt")), 'r') as text_file:
                    for line in text_file:
                        print()
            except:
                report_file.write(f"{line}")


# # Environment settings
# import scanpy as sc
# import os
# import numpy as np
# import pandas as pd
#
# """"  ###### VARIABLES ########  """
#
# h5ad_data_path = '../data/KPTracer-Data_divided/h5ad_files/'
# character_matrix_path = '../data/KPTracer-Data/trees/'
#
# """"  ###### VARIABLES ########  """
#
# # for all interested cell clusters, find the unique tumor classification
# all_files = os.listdir(h5ad_data_path)
# files = [file for file in all_files if not file.startswith('.DS_Store')]
#
# file_path = "report2_refined.txt"
# with open(file_path, 'w') as report_file:
#     for file in files:
#         print("cell cluster: ", file)
#         #report_file.write(f"\ncell cluster {file}")
#
#         final_character_matrix = np.zeros((2,2))
#         lineage_not_found = []
#         adams = sc.read(os.path.join(h5ad_data_path, file))
#         unique_tumor_values = adams.obs['Tumor'].unique()
#         tumor_num =0
#
#         print("num of unique tumor ", len(unique_tumor_values))
#         #report_file.write(f"\nnum of unique tumor {unique_tumor_values}")
#
#         # Yang provided the character matrix as per tumors.
#         # for all tumor values, read the character matrix.
#         for value in unique_tumor_values:
#             tumor_character_matrix = []
#             tumor_file_name = os.path.join(character_matrix_path, (value+"_character_matrix.txt"))
#
#             try: #(some tumor cells do not have lineage data)
#                 with open(os.path.join(tumor_file_name), 'r') as text_file:
#                     # Read the entire contents of the file into a string
#                     for line in text_file:
#                         print()
#             except:
#                 with open(os.path.join(character_matrix_path,
#                 (value.split("_").pop(0)+"_"+value.split("_").pop(1)+"_All_character_matrix.txt")), 'r') as text_file:
#                     for line in text_file:
#                         print()
#                 print(value, "\n")
#                 report_file.write(f"{value }: DONE \n")
#
#             print()