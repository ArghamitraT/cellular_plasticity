
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
def calculate_entropy(probabilities):
    """Calculate entropy for a given probability distribution."""
    #probabilities = np.array([float(p) for p in probabilities])
    probabilities = np.array([p+1e-10 for p in probabilities])
    total_sum = np.sum(probabilities)

    # Calculate percentage
    percentage = 100 * probabilities / total_sum
    entropy = -np.sum(probabilities * np.log(probabilities + 1e-10))  # Avoid log(0) error
    return entropy


def compute_entropy(Gaussian_folder, cell_types_to_consider, query_cell ):
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
        # Likelihoods
        likelihoods = {}
        for cell_type in cell_types_to_consider:
            pdf_value = multivariate_normal.pdf(query_cell,
                        mean=means[cell_type].squeeze(), cov=covariances[cell_type].squeeze())
            likelihoods[cell_type] = pdf_value + epsilon

        # Evidence
        evidence = 0
        for cell_type in cell_types_to_consider:
            evidence += likelihoods[cell_type] * priors[cell_type]

        # Posterior Probabilities
        posterior_probabilities = {}
        for cell_type in cell_types_to_consider:
            posterior_prob = (likelihoods[cell_type] * priors[cell_type]) / evidence
            posterior_probabilities[cell_type] = posterior_prob

        # Convert the dictionary to a long-form DataFrame
        df = pd.DataFrame(posterior_probabilities)
        df_long = df.melt(var_name='Cell Type', value_name='Posterior Probability')

        # Creating the violin plot
        plt.figure(figsize=(12, 6))  # Adjust the size as needed
        sns.violinplot(x='Cell Type', y='Posterior Probability', data=df_long)
        plt.title('Violin Plot of Posterior Probabilities by Cell Type')
        plt.xticks(rotation=45)  # Rotate the x-axis labels for better readability
        plt.show()

        # Entropy
        entropy = []
        post_probability = []
        first_key = cell_types_to_consider[0]
        for i in range(len(posterior_probabilities[first_key])):
            probabilities_at_i = [probabilities[i] for probabilities in posterior_probabilities.values()]
            entropy_at_i = calculate_entropy(probabilities_at_i)
            entropy.append(entropy_at_i)
            post_probability.append(probabilities_at_i)

            print("entropy cell no: ", i)

        return entropy, likelihoods, post_probability


# Term2: Log Density and Distance Multiplication
def calculate_term2(dist_matrix, log_densities, cell_types, query_cell):
    term2_values = []
    for i in range(len(query_cell)):
        weighted_dist = 0
        positive_log_densities = {ct: log_densities[ct][i] for ct in cell_types if log_densities[ct][i] > 0}
        for ct1, ct2 in combinations(positive_log_densities.keys(), 2):
            distance = dist_matrix[cell_types.index(ct1)][cell_types.index(ct2)]
            weighted_dist += distance * positive_log_densities[ct1] * positive_log_densities[ct2]
            # weighted_dist += distance * (positive_log_densities[ct1] + positive_log_densities[ct2])/2
        weighted_dist = weighted_dist
        term2_values.append(weighted_dist)
        print("term2 cell no: ", i)
    return term2_values


def process_clusters(Y_clusters, S_cluster):
    # Load the data
    data_path = '../../data/MSK_tumor_data/adata_combined_mnnc_010920_SCANTN2.h5ad'
    gaussian_path = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/fitted_gaussian_noGMM')
    bhatt_dist_path = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/bhatt_dist_NO_GMM.npy')
    adams_comp = sc.read(data_path)

    # Initialize a dictionary to store results
    cluster_results = {}
    entropy_d = {}
    term2_d = {}
    post_prob_d = {}

    # Iterate over each cluster
    for Y_cluster in Y_clusters:
        # Filter data
        if S_cluster:
            adams1 = adams_comp[adams_comp.obs['Cluster-Name'].isin([Y_cluster])]
            adams1 = adams1[adams1.obs['celltype_hint'].isin([S_cluster])]

        else:
            adams1 = adams_comp[adams_comp.obs['cell_type_med'].isin([Y_cluster])]

        # if S_cluster:
        #     adams1 = adams_comp[adams_comp.obs['Cluster-Name'].isin([Y_cluster])]
        #     adams1 = adams1[adams1.obs['celltype_hint'].isin([S_cluster])]
        #
        # else:
        #     adams1 = adams_comp[adams_comp.obs['Cluster-Name'].isin([Y_cluster])]
        #
        # tumor_names = adams1.obs['Tumor']
        # tumor_names_arr = np.array(tumor_names.unique())
        # KP_tumor_list = []
        # # for tumor in tumor_names_arr:
        # #     try:
        # #         if tumor.split("_")[1] != 'NT':
        # #             KP_tumor_list.append(tumor)
        # #     except:
        # #         print()
        # #
        # # # adams = adams2[adams2.obs.iloc[:, -1].isin(['Mesenchymal-1'])]
        #
        # #(AT) comment
        # KP_tumor_list = ['3433_NT_T2']
        # adams1 = adams1[adams1.obs['Tumor'].isin(KP_tumor_list)]

        query_cell = adams1.obsm['X_scimilarity']

        # Prepare frequency matrix
        freq_mat = utils.make_frequency_matrix(adams1.obs['sc_hits'])
        freq_mat.columns = freq_mat.columns.str.replace(' ', '_')
        cell_types_to_consider = freq_mat.columns.tolist()

        row_sums = np.sum(freq_mat, axis=0)

        # Step 2: Plot the resulting vector in a bar plot
        plt.figure(figsize=(20, 6))  # Adjust the figure size as per your preference
        plt.bar(range(len(row_sums)), row_sums)
        plt.title('Sum of Rows in Frequency Matrix T cell')
        plt.xlabel('Row Index')
        plt.ylabel('Sum')
        plt.show()

        # Compute entropy and likelihoods
        entropy, _, post_prob= compute_entropy(gaussian_path, cell_types_to_consider, query_cell)

        # Calculate log densities
        # log_densities = {cell_type: np.log(likelihoods[cell_type]) for cell_type in cell_types_to_consider}

        # # Calculate Term2
        # dist_matrix = np.load(bhatt_dist_path)
        # term2 = calculate_term2(dist_matrix, log_densities, cell_types_to_consider, query_cell)
        #
        # # Multiply Term1 and Term2
        # final_result = [e * t2 for e, t2 in zip(entropy, term2)]

        term2 = 0
        final_result = 0
        post_prob_d[Y_cluster] = post_prob

        cluster_results[Y_cluster] = final_result
        entropy_d[Y_cluster] = entropy
        term2_d[Y_cluster] = term2

        print("done: ", Y_cluster)

    # return cluster_results, entropy_d, term2_d
    return  post_prob

# Example usage
#Y_clusters = ["AT2-like", "Endoderm-like"]  # Add your cluster names here
#S_cluster = "hepatocyte"

S_cluster = ""
# Y_clusters = ["SCLC-A", "T cell"]
Y_clusters = ["T cell"]
# plasticity_score, entropy_d, term2_d = process_clusters(Y_clusters, S_cluster)
post_prob_d = process_clusters(Y_clusters, S_cluster)
print()

# results = pd.Series(entropy_d)
# plt.figure(figsize=(16, 10))
# colors = ['skyblue', 'lightgreen', 'violet']  # Define colors for each cluster
# sns.violinplot(data=results)
# plt.xticks(ticks=range(len(results.index)), labels=results.index)
# plt.xlabel("Cell clusters")
# plt.ylabel("entropy")
# plt.title("MSK Lung cancer data")
# plt.grid()
# # plt.savefig("../figures/"+utils_AT.create_image_name("entropy"))
# plt.show()

df = pd.DataFrame(post_prob_d)
df_long = df.melt(var_name='Probability Index', value_name='Probability')

# Creating the violin plot
plt.figure(figsize=(20, 10))  # Adjust the figure size as needed
sns.violinplot(x='Probability Index', y='Probability', data=df_long)
plt.title('Violin Plot of Posterior Probabilities SCLC-A')
plt.xticks(rotation=90)  # Rotate the x-axis labels if they are crowded
plt.show()

# results = pd.Series(term2_d)
# plt.figure(figsize=(16, 10))
# colors = ['skyblue', 'lightgreen', 'violet']  # Define colors for each cluster
# sns.violinplot(data=results)
# plt.xticks(ticks=range(len(results.index)), labels=results.index)
# plt.xlabel("Cell clusters")
# plt.ylabel("term 2")
# plt.title("MSK Lung cancer data")
# plt.grid()
# plt.savefig("../figures/"+utils_AT.create_image_name("term2"))
# plt.show()

# results = pd.Series(plasticity_score)
# plt.figure(figsize=(16, 10))
# colors = ['skyblue', 'lightgreen', 'violet']  # Define colors for each cluster
# sns.violinplot(data=results)
# plt.xticks(ticks=range(len(results.index)), labels=results.index)
# plt.xlabel("Cell clusters")
# plt.ylabel("plasticity_score")
# plt.title("MSK Lung cancer data")
# plt.grid()
# plt.savefig("../figures/"+utils_AT.create_image_name("plasticty_score"))
# plt.show()
# print()

# for i, cluster in enumerate(Y_clusters):
#     sns.violinplot(data=results[cluster], color=colors[i])
# Now `results` contains the final results for each cluster in Y_clusters

# Plot for Cluster 1
#plt.subplot(1, 2, 1)  # (number of rows, number of columns, index of the current plot)
#sns.violinplot(data=results, color='skyblue')
#
#sns.violinplot(x = results.index, y=results, color='skyblue')
#plt.xticks(rotation=90)
