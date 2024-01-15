import os
import numpy as np
import pandas as pd
import scanpy as sc
import utils_AT as utils
from scipy.stats import multivariate_normal
import math

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


#def kl_divergence(mu1, sigma1, mu2, sigma2):
    # [The rest of the function remains the same]

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
        return entropy

# Load the data
data_path = '../../data/KPTracer-Data/expression/adata_processed_comp_SCANTN2.h5ad'
adams_comp = sc.read(data_path)

Y_cluster = "Endoderm-like"
S_cluster = "epithelial cell"
# Filter the data
# adams1 = adams_comp[adams_comp.obs['Cluster-Name'].isin(['Endoderm-like'])]
# adams1 = adams1[adams1.obs['celltype_hint'].isin(['secretory cell'])]

adams1 = adams_comp[adams_comp.obs['Cluster-Name'].isin([Y_cluster])]
#S_cluster = adams1.obs['celltype_hint'][0]
adams1 = adams1[adams1.obs['celltype_hint'].isin([S_cluster])]

# # Suppose `query_cell` is the new data point for which you want to calculate the log-likelihood
query_cell = adams1.obsm['X_scimilarity']

# Prepare the frequency matrix
freq_mat = utils.make_frequency_matrix(adams1.obs['sc_hits'])
freq_mat.columns = freq_mat.columns.str.replace(' ', '_')
cell_types_to_consider = freq_mat.columns.tolist()

# Folder paths
folder_path1 = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/fitted_gaussian')
folder_path2 = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/fitted_gaussian_noGMM')

# Compute distances for each folder
entropy1 = pd.Series(compute_entropy(folder_path1, cell_types_to_consider))
entropy2 = pd.Series(compute_entropy(folder_path2, cell_types_to_consider))

utils.violin_plot_figures(entropy1, entropy2,
                    "Fitted_gaussian_W_GMM", "Fitted_gaussian_NO_GMM",
                          ("entropy_"+Y_cluster+"_"+S_cluster))

print()
# If you want to do something with distances_folder1 and distances_folder2, you can do it here.
