import os
import numpy as np
import pandas as pd
import scanpy as sc
import utils_AT as utils
from scipy.stats import multivariate_normal
import math

# Function Definitions

def calculate_entropy(probabilities):
    """Calculate entropy for a given probability distribution."""
    #probabilities = np.array([float(p) for p in probabilities])
    probabilities = np.array([p+1e-10 for p in probabilities])
    total_sum = np.sum(probabilities)

    # Calculate percentage
    percentage = 100 * probabilities / total_sum
    #entropy = -np.sum(probabilities * np.log(probabilities + 1e-10))  # Avoid log(0) error
    entropy = -np.sum(probabilities * np.log(probabilities + 1e-10))  # Avoid log(0) error
    return entropy

# Main Analysis

# Load the data
data_path = '../../data/KPTracer-Data/expression/adata_processed_comp_SCANTN2.h5ad'
adams_comp = sc.read(data_path)

# Filter the data
adams1 = adams_comp[adams_comp.obs['Cluster-Name'].isin(['Endoderm-like'])]
adams1 = adams1[adams1.obs['celltype_hint'].isin(['secretory cell'])]

# Suppose `query_cell` is the new data point for which you want to calculate the log-likelihood
query_cell = adams1.obsm['X_scimilarity']

# Prepare the frequency matrix
freq_mat = utils.make_frequency_matrix(adams1.obs['sc_hits'])
freq_mat.columns = freq_mat.columns.str.replace(' ', '_')
cell_types_to_consider = freq_mat.columns.tolist()

# Gaussian parameters path
Gaussian_folder = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/fitted_gaussian')

# Load means and variances for the cell types
means = {cell_type: np.load(os.path.join(Gaussian_folder, f'{cell_type}_mean.npy'))
         for cell_type in cell_types_to_consider}
covariances = {cell_type: np.load(os.path.join(Gaussian_folder, f'{cell_type}_cov.npy'))
               for cell_type in cell_types_to_consider}

# Read priors
priors = utils.read_csv_column_perRow(
    '../../files/scimilarity_trainingset_cellnum_AuthorLable_original.csv',
    'Cell_Type', 'Number of Samples', cell_types_to_consider)

epsilon = 1e-10  # Small value to prevent log(0)
likelihoods = {cell_type: (multivariate_normal.pdf(query_cell,
            mean=means[cell_type].squeeze(), cov=covariances[cell_type].squeeze()) + epsilon)
               for cell_type in cell_types_to_consider}

evidence = sum(likelihoods[cell_type] * priors[cell_type] for cell_type in cell_types_to_consider)
posterior_probabilities = {cell_type: (likelihoods[cell_type] * priors[cell_type]) / evidence
                           for cell_type in cell_types_to_consider}

# Compute entropy for each data point
entropy = [calculate_entropy([probabilities[i] for probabilities in posterior_probabilities.values()])
           for i in range(len(next(iter(posterior_probabilities.values()))))]

print()





# import matplotlib.pyplot as plt
#
# # Assuming 'cell 1' refers to the first element in the 'query_cell'
# cell_index = 0
#
# # Extract frequencies for cell 1 and convert to percentage
# frequencies_cell1 = freq_mat.iloc[cell_index].values
# total_frequency = np.sum(frequencies_cell1)
# percentages_cell1 = (frequencies_cell1 / total_frequency) * 100
#
# # Extract likelihoods and posterior probabilities for cell 1, and convert to percentage
# likelihoods_cell1 = [likelihoods[cell_type][cell_index] for cell_type in cell_types_to_consider]
# posterior_probs_cell1 = [(posterior_probabilities[cell_type][cell_index] * 100) for cell_type in cell_types_to_consider]
#
# # Sort cell types based on frequency percentages in descending order
# sorted_indices = np.argsort(-percentages_cell1)
# sorted_cell_types = np.array(cell_types_to_consider)[sorted_indices]
#
# # Arrange other metrics according to the sorted order
# sorted_percentages = percentages_cell1[sorted_indices]
# sorted_likelihoods = np.array(likelihoods_cell1)[sorted_indices]
# sorted_posterior_probs = np.array(posterior_probs_cell1)[sorted_indices]
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
#
# # Show plot
# plt.tight_layout()
# plt.savefig(utils.create_image_name('graph_'))
# plt.show()
