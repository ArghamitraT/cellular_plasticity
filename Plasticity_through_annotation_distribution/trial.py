import numpy as np
from itertools import combinations

# ... [Your existing functions and data loading code] ...

def calculate_term2(dist_matrix, likelihoods, cell_types_to_consider, query_cell):
    term2_results = []
    for cell_data in query_cell:  # Iterate over each data point
        # Calculate log density for each cell type
        log_densities = {cell_type: np.log(multivariate_normal.pdf(
            cell_data, mean=means[cell_type].squeeze(), cov=covariances[cell_type].squeeze()) + epsilon)
                         for cell_type in cell_types_to_consider}
        # Filter cell types with positive log density
        positive_log_densities = {ct: log_densities[ct] for ct in cell_types_to_consider if log_densities[ct] > 0}
        # Calculate Term2 for each pair of cell types
        term2_value = sum(dist_matrix[cell_types_to_consider.index(ct1)][cell_types_to_consider.index(ct2)] *
                          positive_log_densities[ct1] * positive_log_densities[ct2]
                          for ct1, ct2 in combinations(positive_log_densities, 2))
        term2_results.append(term2_value)
    return term2_results

# Calculate entropy (Term1)
entropy1 = pd.Series(compute_entropy(folder_path1, cell_types_to_consider))
entropy2 = pd.Series(compute_entropy(folder_path2, cell_types_to_consider))

# Calculate Term2 for each folder
term2_folder1 = calculate_term2(distances_folder1, likelihoods, cell_types_to_consider, query_cell)
term2_folder2 = calculate_term2(distances_folder2, likelihoods, cell_types_to_consider, query_cell)

# Multiply Term1 and Term2
final_result_folder1 = entropy1 * term2_folder1
final_result_folder2 = entropy2 * term2_folder2

# final_result_folder1 and final_result_folder2 contain the product of Term1 and Term2 for each data point in the respective folders


# import matplotlib.pyplot as plt
# import numpy as np
# import os
# import pandas as pd
#
# df = pd.read_csv(os.path.join(os.getcwd(), '../../files/scimilarity_trainingset_cellnum_AuthorLable_original.csv'))
# #cell_types = cell_type_data['Cell_Type']  # Replace 'Cell_Type' with the actual column name
# df['Cell_Type'] = df['Cell_Type'].str.replace(r'.*kidney.*', 'kdny', regex=True)
#
#
# # Example data (replace this with your actual data)
# # data = np.array([0.01, 0.1, 0.3, 0.5, 0.8, 1, 2, 5, 10, 20])
# data = np.linspace(0.1, 20, 4000)  # 400 points between 0.1 and 20
#
#
# # Adding a small constant to avoid log(0)
# epsilon = 1e-10
# data = data + epsilon
#
# # Calculate different logarithmic transformations
# # log2_data = np.log2(data)
# loge_data = np.log(data)  # Natural logarithm
# # log10_data = np.log10(data)
#
# log5_data = np.log(data) / np.log(5)
# log8_data = np.log(data) / np.log(8)  # Natural logarithm
# log10_data = np.log10(data)
#
# # Plotting
# plt.figure(figsize=(10, 6))
#
# # Plot each logarithmic transformation
# plt.plot(data, log5_data,  label='Log5')
# # plt.plot(data, log8_data,  label='Log8')
# plt.plot(data, log10_data,  label='Log10')
# plt.plot(data, loge_data,  label='Loge')
#
# # Adding title, labels, and legend
# plt.title('Logarithmic Transformations')
# plt.xlabel('Original Data')
# plt.ylabel('Transformed Data')
# plt.legend()
# plt.grid()
# plt.savefig('log2.png')
#
# # Show the plot
# plt.show()
# print()
