from scipy.stats import multivariate_normal
import numpy as np
from sklearn.mixture import GaussianMixture
import os
import scanpy as sc


# Suppose `data_point` is the new data point for which you want to calculate the log-likelihood
data_path = '../../data/KPTracer-Data/expression/adata_processed_combined_SCANTN.h5ad'
adams_comp = sc.read(data_path)
adams1 = adams_comp[adams_comp.obs['Cluster-Name'].isin(['Endoderm-like'])]
adams1 = adams1[adams1.obs['celltype_hint'].isin(['secretory cell'])]


data_point = adams1.obsm['X_scimilarity']  # Replace with your actual data point

# Path to the directory with the fitted Gaussian parameters
Gaussian_folder = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/fitted_gaussian')

# Dictionary to store the log-likelihoods for each cell type
log_likelihoods = {}

cell_type = 'secretory_cell'

# Load the mean and covariance
mean = np.load(os.path.join(Gaussian_folder, f'{cell_type}_mean.npy'))
cov = np.load(os.path.join(Gaussian_folder, f'{cell_type}_cov.npy'))
cov = cov.squeeze()

# Create a Gaussian distribution with the loaded parameters
gaussian = multivariate_normal(mean=mean.flatten(), cov=cov)

# Calculate the log-likelihood of the data point under this Gaussian
log_likelihood = gaussian.logpdf(data_point)

# Store the log-likelihood in the dictionary
log_likelihoods[cell_type] = log_likelihood

print()


# Loop over each saved Gaussian parameters file
for file_name in os.listdir(Gaussian_folder):
#for file_name in file_name_arr:
    if file_name.endswith("_mean.npy"):
        cell_type = file_name.replace('_mean.npy', '')

        # Load the mean and covariance
        mean = np.load(os.path.join(Gaussian_folder, f'{cell_type}_mean.npy'))
        cov = np.load(os.path.join(Gaussian_folder, f'{cell_type}_cov.npy'))
        cov = cov.squeeze()

        # Create a Gaussian distribution with the loaded parameters
        gaussian = multivariate_normal(mean=mean.flatten(), cov=cov)

        # Calculate the log-likelihood of the data point under this Gaussian
        log_likelihood = gaussian.pdf(data_point)+1e-250

        # Store the log-likelihood in the dictionary
        log_likelihoods[cell_type] = log_likelihood

        print()
# # Maximum likelihood classification
# max_likelihood = max(likelihoods, key=likelihoods.get)
#
# # Bayesian classification
# bayesian_classification = max(posterior_probabilities, key=posterior_probabilities.get)




# # Dictionary to store the log-likelihoods for each cell type
# likelihoods = {}
#
# # Loop over each saved Gaussian parameters file
# #for file_name in os.listdir(Gaussian_folder):
# for column in freq_mat.columns:
#     cell_type = column.replace(' ', '_')
#     # if file_name.endswith("_mean.npy"):
#     #     cell_type = file_name.replace('_mean.npy', '')
#
#     # Load the mean and covariance
#     mean = np.load(os.path.join(Gaussian_folder, f'{cell_type}_mean.npy'))
#     cov = np.load(os.path.join(Gaussian_folder, f'{cell_type}_cov.npy'))
#     cov = cov.squeeze()
#
#     likelihoods[cell_type] = (multivariate_normal.pdf
#             (query_cell, mean=means[cell_type], cov=covariances[cell_type]))
#
#     # Create a Gaussian distribution with the loaded parameters
#     gaussian = multivariate_normal(mean=mean.flatten(), cov=cov)
#
#     # Calculate the log-likelihood of the data point under this Gaussian
#     log_likelihood = gaussian.pdf(data_point)+1e-250
#
#     # Store the log-likelihood in the dictionary
#     log_likelihoods[cell_type] = log_likelihood
#
#     print()
# Now log_likelihoods will have the log-likelihoods of `data_point` for each fitted Gaussian (cell type)

"""
import numpy as np
from scipy.stats import multivariate_normal

# Mean and covariance for each cell type
means = {
    'epithelial': np.array([...]),
    'endothelial': np.array([...]),
    # ... other cell types
}
covariances = {
    'epithelial': np.array([...]),
    'endothelial': np.array([...]),
    # ... other cell types
}

# Prior probabilities for each cell type
priors = {
    'epithelial': 0.5,  # example prior
    'endothelial': 0.5,
    # ... other cell types
}

# Query cell data
query_cell = np.array([...])  # Your 128-dimensional data here

# Calculate likelihoods
likelihoods = {}
for cell_type in means:
    likelihoods[cell_type] = multivariate_normal.pdf(query_cell, mean=means[cell_type], cov=covariances[cell_type])

# Bayesian approach
posterior_probabilities = {}
evidence = sum(likelihoods[cell_type] * priors[cell_type] for cell_type in likelihoods)
for cell_type in likelihoods:
    posterior_probabilities[cell_type] = (likelihoods[cell_type] * priors[cell_type]) / evidence

# Maximum likelihood classification
max_likelihood = max(likelihoods, key=likelihoods.get)

# Bayesian classification
bayesian_classification = max(posterior_probabilities, key=posterior_probabilities.get)

"""


# cell_type = 'secretory_cell'
#
# # Load the mean and covariance
# mean = np.load(os.path.join(Gaussian_folder, f'{cell_type}_mean.npy'))
# cov = np.load(os.path.join(Gaussian_folder, f'{cell_type}_cov.npy'))
# cov = cov.squeeze()
#
# # Create a Gaussian distribution with the loaded parameters
# gaussian = multivariate_normal(mean=mean.flatten(), cov=cov)
#
# # Calculate the log-likelihood of the data point under this Gaussian
# log_likelihood = gaussian.logpdf(data_point)
#
# # Store the log-likelihood in the dictionary
# log_likelihoods[cell_type] = log_likelihood
#
# print()

"""
from scipy.stats import multivariate_normal
import numpy as np
from sklearn.mixture import GaussianMixture
import os
import scanpy as sc
import utils_AT as utils
import pandas as pd

# Calculate entropy for a given probability distribution
def calculate_entropy(probabilities):
    probabilities = [float(p) for p in probabilities]
    probabilities = np.array(probabilities)
    entropy = -np.sum(probabilities * np.log2(probabilities + 1e-10))  # Adding a small constant to avoid log(0)
    return entropy

# Suppose `data_point` is the new data point for which you want to calculate the log-likelihood
data_path = '../../data/KPTracer-Data/expression/adata_processed_comp_SCANTN2.h5ad'
adams_comp = sc.read(data_path)
adams1 = adams_comp[adams_comp.obs['Cluster-Name'].isin(['Endoderm-like'])]
adams1 = adams1[adams1.obs['celltype_hint'].isin(['secretory cell'])]
query_cell = adams1.obsm['X_scimilarity']  # Replace with your actual data point

freq_mat = utils.make_frequency_matrix(adams1.obs['sc_hits'])
freq_mat.columns = freq_mat.columns.str.replace(' ', '_')
cell_types_to_consider = freq_mat.columns.tolist()

# Path to the directory with the fitted Gaussian parameters
Gaussian_folder = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/fitted_gaussian')
priors = utils.read_csv_column_perRow(
    '../../files/scimilarity_trainingset_cellnum_AuthorLable.csv',
    'Cell_Type', 'Number of Samples', cell_types_to_consider)


# Load means and variances for the selected cell types from numpy files
means = {}
covariances = {}

for cell_type in cell_types_to_consider:
    means[cell_type] = np.load(os.path.join(Gaussian_folder, f'{cell_type}_mean.npy'))
    covariances[cell_type] = np.load(os.path.join(Gaussian_folder, f'{cell_type}_cov.npy'))

# # Prior probabilities for each cell type (you can adjust these as needed)
# priors = {cell_type: 1.0 / len(cell_types_to_consider) for cell_type in cell_types_to_consider}


# Calculate likelihoods
likelihoods = {}
for cell_type in cell_types_to_consider:
    likelihoods[cell_type] = multivariate_normal.pdf(query_cell,
            mean=means[cell_type].squeeze(), cov=covariances[cell_type].squeeze())

# Bayesian approach
posterior_probabilities = {}
evidence = sum(likelihoods[cell_type] * priors[cell_type] for cell_type in cell_types_to_consider)
for cell_type in cell_types_to_consider:
    posterior_probabilities[cell_type] = (likelihoods[cell_type] * priors[cell_type]) / evidence

val = []
_, first_row_values = next(iter(posterior_probabilities.items()))
cell_num = len(first_row_values)

entropy = []
for i in range(cell_num):
    for cell_type, probabilities in posterior_probabilities.items():
        val.append(probabilities[i])
    #val_df = pd.DataFrame({'Value': val})
    entropy.append(calculate_entropy(val))
    i+=1
        #second_values[cell_type] = probabilities[1]  # Extract the second value

print()


"""