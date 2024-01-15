""" This file contains some useful code, used in different reason """
import os
import json
import numpy as np
import pandas as pd
import math

# val1 = math.log(0.9*0.01*0.09)
# #val3 = math.log(0.3*0.3*0.3/3)
# val2 = math.log(0.9*0.01*2*0.02*0.06)
def calculate_entropy(distribution):
    # Ensure no zero probability to avoid log(0)
    #distribution = distribution[distribution > 0]
    return -np.sum(distribution * np.log(distribution))

def calculate_cv(distribution):
    return np.std(distribution) / np.mean(distribution)

def prob_score(dist_q2_lp):
    dist = pd.Series(dist_q2_lp)
    std_q2_lp = dist.std()+1e-200
    #prob_score_q2_lp = (math.prod(dist_q2_lp)*std_q2_lp)
    prob_score_q2_lp = (math.prod(dist_q2_lp))*len(dist_q2_lp)

    return prob_score_q2_lp


# dist_q2_lp = [0.9,0.01,0.01,0.02,0.06]
# prob_score_q2_lp = prob_score(dist_q2_lp)

# dist_q2_mstP = [0.3, 0.2, 0.2, 0.1, 0.2]
# P6 = prob_score(dist_q2_mstP)
# dist_q2_P = [0.6, 0.01, 0.02, 0.01, 0.36]
# P4 = prob_score(dist_q2_P)
# dist_q2_lP = [0.9, 0.01, 0.01, 0.02, 0.06]
# P2 = prob_score(dist_q2_lP)
#
# dist_q1_mstP = [0.33, 0.33, 0.33]
# P5 = prob_score(dist_q1_mstP)
# dist_q1_P = [0.6, 0.3, 0.1]
# P3 = prob_score(dist_q1_P)
# dist_q1_lP = [0.9, 0.01, 0.09]
# P1 = prob_score(dist_q1_lP)

dist_q2_mstP = [0.3, 0.2, 0.2, 0.1, 0.2]
P6 = calculate_entropy(dist_q2_mstP)
dist_q2_P = [0.6, 0.01, 0.02, 0.01, 0.36]
P4 = calculate_entropy(dist_q2_P)
dist_q2_lP = [0.9, 0.01, 0.01, 0.02, 0.06]
P2 = calculate_entropy(dist_q2_lP)

dist_q1_mstP = [0.33, 0.33, 0.33]
P5 = calculate_entropy(dist_q1_mstP)
dist_q1_P = [0.6, 0.3, 0.1]
P3 = calculate_entropy(dist_q1_P)
dist_q1_lP = [0.9, 0.01, 0.09]
P1 = calculate_entropy(dist_q1_lP)

#print(P6, P5, P4, P3, P2, P1)

# dist_q2_mstP = [0.3, 0.2, 0.2, 0.1, 0.2]
# P6 = calculate_cv(dist_q2_mstP)
# dist_q2_P = [0.6, 0.01, 0.02, 0.01, 0.36]
# P4 = calculate_cv(dist_q2_P)
# dist_q2_lP = [0.9, 0.01, 0.01, 0.02, 0.06]
# P2 = calculate_cv(dist_q2_lP)
#
# dist_q1_mstP = [0.33, 0.33, 0.33]
# P5 = calculate_cv(dist_q1_mstP)
# dist_q1_P = [0.6, 0.3, 0.1]
# P3 = calculate_cv(dist_q1_P)
# dist_q1_lP = [0.9, 0.01, 0.09]
# P1 = calculate_cv(dist_q1_lP)

print()

""" for a large csv file print out the row numbers """
# full_path = os.path.join(os.getcwd(), '../../models/query_model_v1/oo_meta.csv')
#
# with open(full_path, 'r') as csvfile:
#     reader = csv.reader(csvfile)
#     for row_number, row in enumerate(reader):
#         print(f'Row {row_number}')
#

""" from the train_meta.csv file find out unique cell types and the indexes for each cell type. Save the dictionary as json """
# full_path = os.path.join(os.getcwd(), '../../models/query_model_v1/train_meta.csv')
# df = pd.read_csv(full_path)
# unique_cell_types = df['author_label'].unique()
#
# # Get the index values for each unique cell type
# index_values_by_type = {cell_type: df.index[df['author_label'] == cell_type].tolist() for cell_type in unique_cell_types}
# for cell_type, indices in index_values_by_type.items():
#     print(f"Indices for {cell_type}: {indices}")
#
# import json
# with open('Index_values_by_type.json', 'w') as file:
#     json.dump(index_values_by_type, file, indent=4)
#
# print()


""" figure MAP, density, % """

# import matplotlib.pyplot as plt
#
# # Assuming 'cell 1' refers to the first element in the 'query_cell'
# cell_index = 0
#
# # Extract frequencies, likelihoods, and posterior probabilities for cell 1
# frequencies_cell1 = freq_mat.iloc[cell_index].values
# likelihoods_cell1 = np.array([likelihoods[cell_type][cell_index] for cell_type in cell_types_to_consider])
# posterior_probs_cell1 = np.array([(posterior_probabilities[cell_type][cell_index] * 100) for cell_type in cell_types_to_consider])
#
# # Convert frequencies to percentage
# total_frequency = np.sum(frequencies_cell1)
# percentages_cell1 = (frequencies_cell1 / total_frequency) * 100
#
# # Filter data where likelihoods are greater than 0
# valid_indices = likelihoods_cell1 > 0
# filtered_cell_types = np.array(cell_types_to_consider)[valid_indices]
# filtered_percentages = percentages_cell1[valid_indices]
# filtered_likelihoods = likelihoods_cell1[valid_indices]
# filtered_posterior_probs = posterior_probs_cell1[valid_indices]
#
# # Sort filtered cell types based on frequency percentages in descending order
# sorted_indices = np.argsort(-filtered_percentages)
# sorted_cell_types = filtered_cell_types[sorted_indices]
# sorted_percentages = filtered_percentages[sorted_indices]
# sorted_likelihoods = filtered_likelihoods[sorted_indices]
# sorted_posterior_probs = filtered_posterior_probs[sorted_indices]
#
# # Prepare data for plotting
# index = np.arange(len(sorted_cell_types))
#
# # Plotting
# plt.figure(figsize=(15, 8))
# plt.scatter(index, sorted_percentages, label='Frequency Percentage', color='r')
# plt.scatter(index, sorted_likelihoods, label='Log Likelihood', color='g')
# plt.scatter(index, sorted_posterior_probs, label='Posterior Probability Percentage', color='b')
#
# # Connect points with lines
# plt.plot(index, sorted_percentages, color='r', alpha=0.5)
# plt.plot(index, sorted_likelihoods, color='g', alpha=0.5)
# plt.plot(index, sorted_posterior_probs, color='b', alpha=0.5)
#
# # Add labels and title
# plt.xlabel('Cell Types (sorted by Frequency Percentage)')
# plt.ylabel('Scores/Percentages')
# plt.title('Frequency Percentages, Likelihood, and Posterior Probabilities (as %) for Cell 1')
# plt.xticks(index, sorted_cell_types, rotation=45)
# plt.legend()
# plt.grid()
#
# # Show plot
# plt.tight_layout()
# plt.savefig(utils.create_image_name('graph_'))
# plt.show()
#
# print()

