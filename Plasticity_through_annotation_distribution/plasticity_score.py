from itertools import combinations
import numpy as np
import pandas as pd
import scanpy as sc
import utils_AT as utils
from scipy.stats import multivariate_normal
import math
import os

# Term1: Calculating Entropy for Each Cell
def calculate_entropy(probabilities):
    """Calculate entropy for a given probability distribution."""
    #probabilities = np.array([float(p) for p in probabilities])
    probabilities = np.array([p+1e-10 for p in probabilities])
    total_sum = np.sum(probabilities)

    # Calculate percentage
    percentage = 100 * probabilities / total_sum
    entropy = -np.sum(probabilities * np.log(probabilities + 1e-10))  # Avoid log(0) error
    return entropy


def compute_entropy(Gaussian_folder, cell_types_to_consider ):
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
                                                           mean=means[cell_type].squeeze(),
                                                           cov=covariances[cell_type].squeeze()) + epsilon)
                       for cell_type in cell_types_to_consider}

        evidence = sum(likelihoods[cell_type] * priors[cell_type] for cell_type in cell_types_to_consider)
        posterior_probabilities = {cell_type: (likelihoods[cell_type] * priors[cell_type]) / evidence
                                   for cell_type in cell_types_to_consider}

        # Compute entropy for each data point
        entropy = [calculate_entropy([probabilities[i] for probabilities in posterior_probabilities.values()])
                   for i in range(len(next(iter(posterior_probabilities.values()))))]
        return entropy, likelihoods


# Term2: Log Density and Distance Multiplication
def calculate_term2(dist_matrix, log_densities, cell_types):
    term2_values = []
    for i in range(len(query_cell)):
        positive_log_densities = {ct: log_densities[ct][i] for ct in cell_types if log_densities[ct][i] > 0}
        for ct1, ct2 in combinations(positive_log_densities.keys(), 2):
            distance = dist_matrix[cell_types.index(ct1)][cell_types.index(ct2)]
            term2_value = distance * positive_log_densities[ct1] * positive_log_densities[ct2]
            term2_values.append(term2_value)
    return term2_values

# paths
data_path = '../../data/KPTracer-Data/expression/adata_processed_comp_SCANTN2.h5ad'
gaussian_path = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/fitted_gaussian_noGMM')
bhatt_dist_path = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/bhatt_dist_NO_GMM.npy')

# Load the data
adams_comp = sc.read(data_path)
Y_cluster = "Endoderm-like"
S_cluster = "epithelial cell"
adams1 = adams_comp[adams_comp.obs['Cluster-Name'].isin([Y_cluster])]
#S_cluster = adams1.obs['celltype_hint'][0]
adams1 = adams1[adams1.obs['celltype_hint'].isin([S_cluster])]
query_cell = adams1.obsm['X_scimilarity']

# Prepare the frequency matrix
freq_mat = utils.make_frequency_matrix(adams1.obs['sc_hits'])
freq_mat.columns = freq_mat.columns.str.replace(' ', '_')
cell_types_to_consider = freq_mat.columns.tolist()

entropy, likelihoods = pd.Series(compute_entropy(gaussian_path, cell_types_to_consider))


# Calculate log densities
log_densities = {cell_type: np.log(likelihoods[cell_type]) for cell_type in cell_types_to_consider}

# Calculate term2
dist_matrix = np.load(bhatt_dist_path)
term2 = calculate_term2(dist_matrix, log_densities, cell_types_to_consider)

# Final Step: Multiply Term1 and Term2
final_result = [e * t2 for e, t2 in zip(entropy, term2)]

# final_result contains the product of Term1 and Term2 for each cell


