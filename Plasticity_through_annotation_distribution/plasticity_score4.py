
import os
import numpy as np
import pandas as pd
import scanpy as sc

import utils_AT
import utils_AT as utils
from scipy.stats import multivariate_normal
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns

# Term1: Calculating Entropy for Each Cell
def calculate_entropy(cell_probabilities):
    """Calculate entropy for a given probability distribution."""

    entropies = []

    for cell in cell_probabilities:
        # Ensure no zero probabilities and convert to numpy array
        probabilities = np.array([p + 1e-10 for p in cell])

        # Calculate entropy
        entropy = -np.sum(probabilities * np.log(probabilities))
        entropies.append(entropy)

    return entropies


# Term2: Log Density and Distance Multiplication
def calculate_term2(dist_matrix, log_densities, cell_types, query_cell):

    term2_values = []
    for i in range(len(query_cell)):
        weighted_dist = 0
        positive_log_densities = {ct: log_densities[ct][i] for ct in cell_types if log_densities[ct][i] > 0}
        for ct1, ct2 in combinations(positive_log_densities.keys(), 2):
            distance = dist_matrix[cell_types.index(ct1)][cell_types.index(ct2)]
            weighted_dist += distance * positive_log_densities[ct1] * positive_log_densities[ct2]

        term2_values.append(weighted_dist)
        print("term2 cell no: ", i)
    return term2_values


def process_clusters(Y_clusters, S_cluster):
    data_path = '../../data/MSK_tumor_data/adata_combined_mnnc_010920_SCANTN2.h5ad'
    gaussian_path = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/fitted_gaussian_noGMM')
    bhatt_dist_path = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/bhatt_dist_NO_GMM.npy')
    adams_comp = sc.read(data_path)
    epsilon = 1e-10

    # Initialize a dictionary to store results
    plasticity_score_d = {}
    entropy_d = {}
    term2_d = {}

    for Y_cluster in Y_clusters:

        adams1 = adams_comp[adams_comp.obs['cell_type_med'].isin([Y_cluster])][:5]
        freq_mat = utils.make_frequency_matrix(adams1.obs['sc_hits'])
        freq_mat.columns = freq_mat.columns.str.replace(' ', '_')
        cell_types_to_consider = freq_mat.columns.tolist()

        # load the fitted gaussian mean, variance an cell numbers
        means = {cell_type: np.load(os.path.join(gaussian_path, f'{cell_type}_mean.npy'))
                 for cell_type in cell_types_to_consider}
        covariances = {cell_type: np.load(os.path.join(gaussian_path, f'{cell_type}_cov.npy'))
                       for cell_type in cell_types_to_consider}
        priors = utils.read_csv_column_perRow(
            '../../files/scimilarity_trainingset_cellnum_AuthorLable_original.csv',
            'Cell_Type', 'Number of Samples', cell_types_to_consider)

        final_results_per_cluster = []  # Store final results for each cell in this cluster

        cell = adams1.obsm['X_scimilarity']

        # 2. Calculate Likelihood
        likelihoods = {}
        for cell_type in cell_types_to_consider:
            pdf_value = multivariate_normal.pdf(cell,
                                                mean=means[cell_type].squeeze(),
                                                cov=covariances[cell_type].squeeze())
            likelihoods[cell_type] = pdf_value + epsilon

        # Evidence for this cell
        evidence = sum(likelihoods[cell_type] * priors[cell_type] for cell_type in cell_types_to_consider)


        # 4. Calculate Posterior Probability
        posterior_probabilities = {ct: (likelihoods[ct] * priors[ct]) / evidence for ct in
                                   cell_types_to_consider}

        # fetch information for every cell
        num_cells = len(next(iter(posterior_probabilities.values())))
        # Initialize a dictionary to store cell-wise data
        cell_wise_data = [[] for _ in range(num_cells)]

        # Iterate over each row and distribute its data to the respective cells
        for row in posterior_probabilities.values():
            for cell_index, cell_value in enumerate(row):
                cell_wise_data[cell_index].append(cell_value)

        # 5. Calculate Entropy (Term 1)
        cell_wise_prob = [np.array(cell) for cell in cell_wise_data]
        entropy_arr = calculate_entropy(cell_wise_prob)

        # 6. Calculate Term 2
        log_densities = {ct: np.log(likelihoods[ct] + epsilon) for ct in cell_types_to_consider}
        dist_matrix = np.load(bhatt_dist_path)
        term2_arr = calculate_term2(dist_matrix, log_densities, cell_types_to_consider, cell)

        # 7. Calculate Final Result
        plasticity_score_arr = [e * t2 for e, t2 in zip(entropy_arr, term2_arr)]

        # create dictionary for every Y cluster
        entropy_d[Y_cluster] = entropy_arr
        term2_d[Y_cluster] = term2_arr
        plasticity_score_d[Y_cluster] = plasticity_score_arr

    return entropy_d, term2_d, plasticity_score_d


S_cluster = ""
Y_clusters = ["SCLC-A", "T cell"]
#plasticity_score_d, entropy_d, term2_d = process_clusters(Y_clusters, S_cluster)
entropy_d, term2_d, plasticity_score_d = process_clusters(Y_clusters, S_cluster)
print()

results = pd.Series(term2_d)
plt.figure(figsize=(16, 10))
colors = ['skyblue', 'lightgreen', 'violet']  # Define colors for each cluster
sns.violinplot(data=results)
plt.xticks(ticks=range(len(results.index)), labels=results.index)
plt.xlabel("Cell clusters")
plt.ylabel("term 2")
plt.title("MSK Lung cancer data")
plt.grid()
plt.savefig("../figures/"+utils_AT.create_image_name("term2"))
plt.show()

results = pd.Series(entropy_d)
plt.figure(figsize=(16, 10))
colors = ['skyblue', 'lightgreen', 'violet']  # Define colors for each cluster
sns.violinplot(data=results)
plt.xticks(ticks=range(len(results.index)), labels=results.index)
plt.xlabel("Cell clusters")
plt.ylabel("entropy")
plt.title("MSK Lung cancer data")
plt.grid()
plt.savefig("../figures/"+utils_AT.create_image_name("entropy"))
plt.show()

results = pd.Series(plasticity_score_d)
plt.figure(figsize=(16, 10))
colors = ['skyblue', 'lightgreen', 'violet']  # Define colors for each cluster
sns.violinplot(data=results)
plt.xticks(ticks=range(len(results.index)), labels=results.index)
plt.xlabel("Cell clusters")
plt.ylabel("plasticity_score")
plt.title("MSK Lung cancer data")
plt.grid()
plt.savefig("../figures/"+utils_AT.create_image_name("plasticty_score"))
plt.show()
print()

# for i, cluster in enumerate(Y_clusters):
#     sns.violinplot(data=results[cluster], color=colors[i])
# Now `results` contains the final results for each cluster in Y_clusters

# Plot for Cluster 1
#plt.subplot(1, 2, 1)  # (number of rows, number of columns, index of the current plot)
#sns.violinplot(data=results, color='skyblue')
#
#sns.violinplot(x = results.index, y=results, color='skyblue')
#plt.xticks(rotation=90)
