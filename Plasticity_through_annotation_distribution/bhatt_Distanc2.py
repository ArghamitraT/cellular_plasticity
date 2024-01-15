import numpy as np
import os
import pandas as pd
import utils_AT
import scipy.cluster.hierarchy as sch
import seaborn as sns
import matplotlib.pyplot as plt

def plot_heatmap_with_dendrogram(dist_matrix, labels, title):
    # Perform hierarchical clustering
    linkage = sch.linkage(sch.distance.squareform(dist_matrix), method='average')

    # Create a dendrogram
    dendro = sch.dendrogram(linkage, labels=labels, leaf_rotation=90)

    # Reorder the distances matrix and labels according to the dendrogram
    ordered_labels = [labels[i] for i in dendro['leaves']]
    ordered_matrix = dist_matrix[:, dendro['leaves']][dendro['leaves'], :]

    # Plot the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(ordered_matrix, cmap='viridis', xticklabels=ordered_labels, yticklabels=ordered_labels)
    plt.title(title)
    plt.show()

def bhattacharyya_distance(mean1, cov1, mean2, cov2):
    # Calculate the average covariance matrix
    Sigma = 0.5 * (cov1 + cov2)

    mean1 = mean1.squeeze()
    mean2 = mean2.squeeze()

    # Compute the Bhattacharyya distance
    diff_mean = mean2 - mean1
    term1 = 0.125 * np.dot(np.dot(diff_mean.T, np.linalg.inv(Sigma)), diff_mean)
    # term2 = 0.5 * np.log(np.linalg.det(Sigma) / np.sqrt(np.linalg.det(cov1) * np.linalg.det(cov2)))
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

    print("running")
    return term1 + term2


#def kl_divergence(mu1, sigma1, mu2, sigma2):
    # [The rest of the function remains the same]

def compute_distances(folder_path, cell_types):
    # Initialize a matrix to store the distances
    num_cell_types = len(cell_types)
    distances = np.zeros((num_cell_types, num_cell_types))

    # Compute pairwise distances
    for i, type1 in enumerate(cell_types):
        mean1 = np.load(os.path.join(folder_path, f'{type1}_mean.npy')).squeeze()
        cov1 = np.load(os.path.join(folder_path, f'{type1}_cov.npy')).squeeze()

        for j, type2 in enumerate(cell_types):
            if j > i:  # To avoid redundant computations
                mean2 = np.load(os.path.join(folder_path, f'{type2}_mean.npy')).squeeze()
                cov2 = np.load(os.path.join(folder_path, f'{type2}_cov.npy')).squeeze()

                distance = bhattacharyya_distance(mean1, cov1, mean2, cov2)
                # Or use KL divergence if preferred
                # distance = kl_divergence(mean1, cov1, mean2, cov2)

                distances[i, j] = distance
                distances[j, i] = distance  # Since the distance is symmetric

    return distances

# Folder paths
folder_path1 = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/fitted_gaussian')
folder_path2 = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/fitted_gaussian_noGMM')

cell_type_data = pd.read_csv(os.path.join(os.getcwd(), '../../files/scimilarity_trainingset_cellnum_AuthorLable_original.csv'))
cell_types = cell_type_data['Cell_Type'].tolist()  # Replace 'Cell_Type' with the actual column name

# Compute distances for each folder
distances_folder1 = compute_distances(folder_path1, cell_types)
distances_folder2 = compute_distances(folder_path2, cell_types)
