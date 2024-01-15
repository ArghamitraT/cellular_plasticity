"""
The code goes to Yang character matrices separater by tumors and
form a character matrix file with all cells
output: save the character matrix as a pkl and csv file
"""

# Environment settings
import os
import numpy as np
import pandas as pd

character_matrix_path = '../data/KPTracer-Data/trees/'

# for all interested cell clusters, find the unique tumor classification
all_files = os.listdir(character_matrix_path)
files = [file for file in all_files if not file.startswith('.DS_Store')]

final_character_matrix = np.zeros((2,2))
tumor_num = 0

# for all files in the tree folder; these files are separated by tumor
for file in files:

    # there are .pkl and .nwk files in the folder, but we only need .txt files
    extension = file.split(".")[-1]
    if extension=='txt':
        tumor_file_name = os.path.join(character_matrix_path, file)

        with open(tumor_file_name, 'r') as text_file:
            tumor_character_matrix = []
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
                line_num +=1

            # in the tumor character matrix, put the tumor name as well
            try:
                tumor_name_separated = tumor_file_name.split('/')[-1].split('.')[0].split('_')
                tumor_name = tumor_name_separated[0]+'_'+tumor_name_separated[1]+'_'+tumor_name_separated[2]
            except:
                print()
            tumor_character_matrix = np.array(tumor_character_matrix)
            tumor_column = np.full((tumor_character_matrix.shape[0], 1), tumor_name)
            tumor_character_matrix = np.concatenate((tumor_column, tumor_character_matrix), axis=1)

            final_cell_clm_num = np.shape(final_character_matrix)[1]
            current_cell_clm_num = np.shape(tumor_character_matrix)[1]

            # merge the tumor character matrix in the final file with all character matrix
            if tumor_num > 0:
                if (final_cell_clm_num > current_cell_clm_num):  # if previous tumor cells have more mutation state
                    matrix_padding = np.full((np.shape(tumor_character_matrix)[0], abs(final_cell_clm_num - current_cell_clm_num)), 0)
                    tumor_character_matrix = np.hstack((tumor_character_matrix, matrix_padding))
                elif (final_cell_clm_num < current_cell_clm_num):  # if current tumor cells have more mutation state
                    matrix_padding = np.full((np.shape(final_character_matrix)[0], abs(final_cell_clm_num - current_cell_clm_num)), 0)
                    final_character_matrix = np.hstack((final_character_matrix, matrix_padding))

                final_character_matrix = np.concatenate((final_character_matrix, tumor_character_matrix),
                                                        axis=0)
            else:  # for first tumor appending
                final_character_matrix = tumor_character_matrix

            print(f"file number {tumor_num}")
            tumor_num += 1

# panda formatting, change the cells as row header
final_character_matrix = pd.DataFrame(final_character_matrix)
header_row = final_character_matrix.iloc[:,1]
#final_character_matrix = final_character_matrix.drop(1, axis=1)
final_character_matrix.set_index(header_row.values, inplace=True)

# save
final_name = "../data/KPTracer-Data_divided/character_matrix/ALL_CELL_character_matrix_2.pkl"
final_character_matrix.to_pickle(final_name)
csv_file_path = "../data/KPTracer-Data_divided/character_matrix/ALL_CELL_character_matrix_2.csv"
final_character_matrix.to_csv(csv_file_path, index=False)

#final_character_matrix.to_pickle('../data/KPTracer-Data_divided/h5ad_filesWLineage/yang_PrEMT_scrtrCell_character_matrix.pkl')