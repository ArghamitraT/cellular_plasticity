import numpy as np
import os
import pandas as pd

def bhattacharyya_distance(mean1, cov1, mean2, cov2):
    # Calculate the average covariance matrix
    Sigma = 0.5 * (cov1 + cov2)

    mean1 = mean1.squeeze()
    mean2 = mean2.squeeze()

    # Compute the Bhattacharyya distance
    diff_mean = mean2 - mean1
    term1 = 0.125 * np.dot(np.dot(diff_mean.T, np.linalg.inv(Sigma)), diff_mean)
    #term2 = 0.5 * np.log(np.linalg.det(Sigma) / np.sqrt(np.linalg.det(cov1) * np.linalg.det(cov2)))
    sign1, logdet1 = np.linalg.slogdet(cov1)
    sign2, logdet2 = np.linalg.slogdet(cov2)
    sign, logdet = np.linalg.slogdet(Sigma)

    # Make sure that all determinants are positive (sign should be 1)
    if sign1 == sign2 == sign == 1:
        term2 = 0.5 * (logdet - 0.5 * (logdet1 + logdet2))
    else:
        # Handle the case where the determinant is negative or zero
        # This might indicate an issue with the covariance matrices
        print("Error: Non-positive determinant encountered.")

    return term1 + term2


def kl_divergence(mu1, sigma1, mu2, sigma2):
    """
    Compute the KL divergence between two multivariate Gaussian distributions.

    Parameters:
    mu1 : array-like
        Mean of the first Gaussian distribution.
    sigma1 : array-like
        Covariance matrix of the first Gaussian distribution.
    mu2 : array-like
        Mean of the second Gaussian distribution.
    sigma2 : array-like
        Covariance matrix of the second Gaussian distribution.

    Returns:
    float
        The KL divergence between the two distributions.
    """
    k = len(mu1)
    sigma2_inv = np.linalg.inv(sigma2)
    trace_term = np.trace(sigma2_inv @ sigma1)
    mean_diff = mu2 - mu1
    mean_term = mean_diff.T @ sigma2_inv @ mean_diff
    det_term = np.log(np.linalg.det(sigma2) / np.linalg.det(sigma1))

    return 0.5 * (trace_term + mean_term - k + det_term)


folder_path = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/fitted_gaussian')
cell_types = cell_type_data['Cell_Type'].tolist()  # Replace 'cell_type' with the actual column name

# Initialize a matrix to store the distances
num_cell_types = len(cell_types)
distances = np.zeros((num_cell_types, num_cell_types))

# Compute pairwise distances
for i, type1 in enumerate(cell_types):
    mean1 = np.load(os.path.join(folder_path, f'{type1}_mean.npy')).squeeze()
    cov1 = np.load(os.path.join(folder_path, f'{type1}_cov.npy'))
    cov1 = cov1.squeeze()

    for j, type2 in enumerate(cell_types):
        if j > i:  # To avoid redundant computations
            mean2 = np.load(os.path.join(folder_path, f'{type2}_mean.npy')).squeeze()
            cov2 = np.load(os.path.join(folder_path, f'{type2}_cov.npy'))
            cov2 = cov2.squeeze()

            distance = bhattacharyya_distance(mean1, cov1, mean2, cov2)

            #distance = kl_divergence(mean1, cov1, mean2, cov2)
            distances[i, j] = distance
            distances[j, i] = distance  # Since the distance is symmetric

        print()
# Example usage
mean1 = np.array([1, 2])
cov1 = np.array([[1, 0], [0, 1]])

mean2 = np.array([4, 5])
cov2 = np.array([[2, 0], [0, 2]])

distance = bhattacharyya_distance(mean1, cov1, mean2, cov2)
print("Bhattacharyya Distance:", distance)
