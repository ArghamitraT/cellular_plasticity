""" annoation distribution, Edited the code to handle three cluster side by side"""

""" environment """
import scanpy as sc
import warnings
import utils_AT
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import json

""" variables """

# need to be Yang's cluster
cluster1 = 'Mesenchymal-1'
cluster2 = 'Pre-EMT'
cluster3 = 'Early gastric'

# cluster1 = 'Mesenchymal-2'
# cluster2 = 'AT1-like'
# cluster3 = 'Early gastric'


# need to be scimilarity cluster
common_cluster = 'NA'
# # need to be Yang's cluster
# cluster1 = 'Mesenchymal-1 (Met)'
# cluster2 = 'Mesenchymal-1'
#
# # need to be scimilarity cluster
# common_cluster = 'fibroblast'

# # need to be Yang's cluster
# cluster1 = 'High plasticity'
# cluster2 = 'AT1-like'
#
# # need to be scimilarity cluster
# common_cluster = 'epithelial cell'

# # need to be Yang's cluster
# cluster1 = 'High plasticity'
# cluster2 = 'Endoderm-like'
#
# # need to be scimilarity cluster
# common_cluster = 'epithelial cell'
title_avg_antn = 'Distribution of Average Percentage '
title_STD = 'Distribution of STD '

""" input-> distance matrix, output-> sigmoid distance matrix """
def adjusted_sigmoid(matrix, k=1):
    return 1 / (1 + np.exp(-k * matrix))


""" average distance multiplied by cell type annotation fraction """
def weighted_annotation(avg_dist, freq):
    # avg distance per celltype is multiplied by percentage of cell type. so Weighted distance for every cell type
    wghtd_antn = avg_dist * (freq / 50)

    # for every score adding the weighted distance
    wghtd_score = wghtd_antn.sum(axis=1)
    return wghtd_score


""" input-> series data of cells and cluster name, output-> frequency matrix of the cell cluster """
# def make_frequency_matrix(celltype_hits_series, dist=0):
#     # Initialize an empty DataFrame
#     frequency_df = pd.DataFrame()
#
#     # Iterate through each row and update the DataFrame
#     for idx, celltype_str in celltype_hits_series.items():
#         celltype_dict = json.loads(celltype_str)
#
#         for cell, freq in celltype_dict.items():
#             if cell not in frequency_df.columns:
#                 frequency_df[cell] = 0
#             if dist:
#                 frequency_df.at[idx, cell] = (freq)
#             else:
#                 frequency_df.at[idx, cell] = int(freq)
#
#     # Fill NaN values with 0
#     frequency_df = frequency_df.fillna(0)
#     return frequency_df
def make_frequency_matrix(celltype_hits_series, dist=0):
    data_dict = {}  # Dictionary to store all data

    # Iterate through each row and populate the data_dict
    for index, (idx, celltype_str) in enumerate(celltype_hits_series.items()):
        celltype_dict = json.loads(celltype_str)

        for cell, freq in celltype_dict.items():
            if cell not in data_dict:
                data_dict[cell] = [0] * len(celltype_hits_series)

            if dist:
                data_dict[cell][index] = (freq)
            else:
                data_dict[cell][index] = int(freq)

    # Convert the data_dict to a DataFrame
    frequency_df = pd.DataFrame(data_dict, index=celltype_hits_series.index)

    # This might not be necessary now, but just in case
    frequency_df = frequency_df.fillna(0)

    return frequency_df


""" input-> weighted annotaion, output-> entropy of the matrix """
def compute_entropy_mat(df):
    # Ensure the DataFrame is in the desired format
    if not isinstance(df, pd.DataFrame):
        raise ValueError("Input should be a pandas DataFrame")

    # Normalize each row
    row_sums = df.sum(axis=1)
    normalized_matrix = df.div(row_sums, axis=0)

    # Compute entropy for each row
    entropies = -np.sum(normalized_matrix * np.log(normalized_matrix + 1e-10),
                        axis=1)  # Added small value to avoid log(0)
    return entropies


""" average distance for every cell type """
# for similar cell type we have added the distances. We need to take the average. For example if we have 5 epithelial cells in to nearest neighbor, all 5 distances are added here. So distance need to be divided by 5 to get the average distance for epithelial cell (discreetizing 1)
def make_avg_dist(frequency_mat, dist_mat, scaling_factor=1):
    x = dist_mat*scaling_factor
    x = x/frequency_mat
    x.fillna(0, inplace=True)
    return x

""" function: input-> frequency matrix, output-> statistical output: avg, std, %"""
def make_summary_df(frequency_df):
    column_averages = frequency_df.mean()
    column_std = frequency_df.std()
    summary_df = pd.DataFrame({
        'Average': column_averages,
        'STD': column_std
    })
    summary_df['Average Percentage (%)'] = (summary_df['Average'] / summary_df['Average'].sum()) * 100
    sorted_summary_df = summary_df.sort_values(by=['Average Percentage (%)'], ascending=False)
    return sorted_summary_df


""" function: input-> normalized annotation distribution, output-> entropy of the distribution  """
def entropy(probabilities):
    """
    Compute the entropy of a normalized distribution.

    Parameters:
    - probabilities: list or numpy array of probabilities. Should sum to 1.

    Returns:
    - Entropy of the distribution.
    """
    # Ensure the probabilities are a numpy array
    probabilities = np.array(probabilities)

    # Filter out zero probabilities to avoid log(0)
    non_zero_probs = probabilities[probabilities > 0]

    # Compute entropy
    H = -np.sum(non_zero_probs * np.log(non_zero_probs))

    return H


import utils_AT as util


def plot_figures(cluster1_data, cluster2_data, cluster3_data, y_lable):
    plt.figure(figsize=(16, 6))

    # Plot for Cluster 1
    plt.subplot(1, 3, 1)  # (number of rows, number of columns, index of the current plot)
    sns.violinplot(data=cluster1_data, color='skyblue')
    plt.xticks(rotation=90)
    plt.xlabel(cluster1 + " " + common_cluster)
    plt.ylabel(y_lable)
    plt.grid()
    plt.title("Violin plot for Cluster 1 (only tumor cells)")

    # Plot for Cluster 2
    plt.subplot(1, 3, 2)
    sns.violinplot(data=cluster2_data, color='lightgreen')
    plt.xticks(rotation=90)
    plt.xlabel(cluster2 + " " + common_cluster)
    plt.ylabel(y_lable)
    plt.grid()
    plt.title("Violin plot for Cluster 2 (only tumor cells)")

    # Plot for Cluster 3
    plt.subplot(1, 3, 3)
    sns.violinplot(data=cluster3_data, color='yellow')
    plt.xticks(rotation=90)
    plt.xlabel(cluster3 + " " + common_cluster)
    plt.ylabel(y_lable)
    plt.grid()
    plt.title("Violin plot for Cluster 3 (only tumor cells)")

    plt.tight_layout()  # Ensures that the plots do not overlap
    plt.savefig(('figures/' + util.create_image_name(y_lable, format='.jpg')), bbox_inches='tight', dpi=150)
    plt.show()

""" read data """
data_path = '../data/KPTracer-Data/expression/adata_processed_comp_SCANTN2.h5ad' #(AT)
adams_comp = sc.read(data_path)
print()

""" Just take the tumor cells """
tumor_names = adams_comp.obs['Tumor']
tumor_names_arr = np.array(tumor_names.unique())
KP_tumor_list = []
for tumor in tumor_names_arr:
     try:
          if tumor.split("_")[1] != 'NT':
               KP_tumor_list.append(tumor)
     except:
          print()

# adams = adams2[adams2.obs.iloc[:, -1].isin(['Mesenchymal-1'])]
adams_comp_original = adams_comp
adams_comp = adams_comp[adams_comp.obs.iloc[:, 6].isin(KP_tumor_list)]

""" takes two cell clusters and make frequency matrix """
# cluster 1
adams1 = adams_comp[adams_comp.obs['Cluster-Name'].isin([cluster1])]

# comment this line

adams11 = adams1
celltype_hits_series1 = adams11.obs['sc_hits']
frequency_df1 = make_frequency_matrix(celltype_hits_series1)
celltype_dist_series1 = adams11.obs['sc_dist_percelltype']
dist_df1 = make_frequency_matrix(celltype_dist_series1, dist=1)

# cluster 2
adams2 = adams_comp[adams_comp.obs['Cluster-Name'].isin([cluster2])]
adams22 = adams2
#celltype_hits_series2 = adams22.obs['celltype_hits']
celltype_hits_series2 = adams22.obs['sc_hits']
frequency_df2 = make_frequency_matrix(celltype_hits_series2)
celltype_dist_series2 = adams22.obs['sc_dist_percelltype']
dist_df2 = make_frequency_matrix(celltype_dist_series2, dist=1)

# cluster 3
adams3 = adams_comp[adams_comp.obs['Cluster-Name'].isin([cluster3])]
adams33 = adams3
#celltype_hits_series2 = adams22.obs['celltype_hits']
celltype_hits_series3 = adams33.obs['sc_hits']
frequency_df3 = make_frequency_matrix(celltype_hits_series3)
celltype_dist_series3 = adams33.obs['sc_dist_percelltype']
dist_df3 = make_frequency_matrix(celltype_dist_series3, dist=1)

print()

""" make statistics """
summary_df1 = make_summary_df(frequency_df1)
summary_df2 = make_summary_df(frequency_df2)
summary_df3 = make_summary_df(frequency_df3)
print()

# """ Approach 2: entropy of the distribution """
# print()
# entropy_clstr1 = entropy(summary_df1['Average Percentage (%)']/100)
# entropy_clstr2 = entropy(summary_df2['Average Percentage (%)']/100)
# entropy_clstr3 = entropy(summary_df3['Average Percentage (%)']/100)
# print("Approach 2: distribution entropy :", cluster1, " ", common_cluster, entropy_clstr1)
# print("Approach 2: distribution entropy :", cluster2, " ", common_cluster, entropy_clstr2)
# print("Approach 2: distribution entropy :", cluster3, " ", common_cluster, entropy_clstr3)

""" Visualize (Approach2) one way"""
cell_entropy_clstr1 = compute_entropy_mat(frequency_df1)
cell_entropy_clstr2 = compute_entropy_mat(frequency_df2)
cell_entropy_clstr3 = compute_entropy_mat(frequency_df3)

""" making a box plot similar to yang """
df1 = cell_entropy_clstr1.reset_index(drop=True).to_frame(name='avg(Entropy_per_cell)')
df1['Cell Type'] = cluster1
df2 = cell_entropy_clstr2.reset_index(drop=True).to_frame(name='avg(Entropy_per_cell)')
df2['Cell Type'] = cluster2
df3 = cell_entropy_clstr3.reset_index(drop=True).to_frame(name='avg(Entropy_per_cell)')
df3['Cell Type'] = cluster3
final_df = pd.concat([df1, df2, df3], ignore_index=True)
plt.figure(figsize=(10, 6))
sns.boxplot(y='Cell Type', x='avg(Entropy_per_cell)', data=final_df)
plt.xlabel('avg(Entropy_per_cell)')
plt.ylabel('Cell Type')
plt.title('Boxplot of avg(Entropy_per_cell) by Cell Type (only tumor cells)')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.tight_layout()
#plt.savefig(('figures/' + util.create_image_name("box_plot", format='.jpg')), bbox_inches='tight', dpi=150)
plt.show()

#plot_figures(cell_entropy_clstr1, cell_entropy_clstr2, cell_entropy_clstr3, y_lable='Entropy_per_cell')

print("Approach 2b: mean entropy :", cluster1, " ", common_cluster, round(cell_entropy_clstr1.mean(), 3))
print("Approach 2b: mean entropy :", cluster2, " ", common_cluster, round(cell_entropy_clstr2.mean(), 3))
print("Approach 2b: mean entropy :", cluster3, " ", common_cluster, round(cell_entropy_clstr3.mean(), 3))

""" Approach 3: weighted annotation """
# for similar cell type we have added the distances. We need to take the average. For example if we have 5 epithelial cells in to nearest neighbor, all 5 distances are added here. So distance need to be divided by 5 to get the average distance for epithelial cell (discreetizing 1)
dist_df11 = make_avg_dist(frequency_df1, dist_df1)
dist_df22 = make_avg_dist(frequency_df2, dist_df2)
dist_df33 = make_avg_dist(frequency_df3, dist_df3)

clstr1_wghtd_score = weighted_annotation(dist_df11, frequency_df1)
clstr2_wghtd_score = weighted_annotation(dist_df22, frequency_df2)
clstr3_wghtd_score = weighted_annotation(dist_df33, frequency_df3)

""" making a box plot similar to yang """
df1 = clstr1_wghtd_score.reset_index(drop=True).to_frame(name='avg(Entropy_per_cell)')
df1['Cell Type'] = cluster1
df2 = clstr2_wghtd_score.reset_index(drop=True).to_frame(name='avg(Entropy_per_cell)')
df2['Cell Type'] = cluster2
df3 = clstr3_wghtd_score.reset_index(drop=True).to_frame(name='avg(Entropy_per_cell)')
df3['Cell Type'] = cluster3
final_df = pd.concat([df1, df2, df3], ignore_index=True)
plt.figure(figsize=(10, 6))
sns.boxplot(y='Cell Type', x='avg(Entropy_per_cell)', data=final_df)
plt.xlabel('avg(weighted_annotation_score_per_cell)')
plt.ylabel('Cell Type')
plt.title('Boxplot of avg(weighted_annotation_score_per_cell) by Cell Type (only tumor cells)')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.savefig(('figures/' + util.create_image_name("box_plot", format='.jpg')), bbox_inches='tight', dpi=150)
plt.show()

""" Visualize (Approach 3) """
plot_figures(clstr1_wghtd_score, clstr2_wghtd_score, clstr3_wghtd_score, y_lable='weighted_annotation_score_per_cell')

# the weighted distance for a cluster
print("Approach 3: weighted annotation :", cluster1, " ", common_cluster, round(clstr1_wghtd_score.mean(),3))
print("Approach 3: weighted annotation :", cluster2, " ", common_cluster, round(clstr2_wghtd_score.mean(),3))
print("Approach 3: weighted annotation :", cluster3, " ", common_cluster, round(clstr3_wghtd_score.mean(),3))
print()
#
# """ Approach 4: sigmoid weighted annotation """
# # k_mat = [1, 0.8, 0.6, 0.4, 0.2]
# # k_mat = [0.1, 0.08, 0.06]
# k_mat = [5]
#
# for k_val in k_mat:
#     # get sigmoid of average distance
#     sig_dist1 = adjusted_sigmoid(dist_df11, k=k_val)
#     sig_dist2 = adjusted_sigmoid(dist_df22, k=k_val)
#
#     sig_clstr1_wghtd_score = weighted_annotation(sig_dist1, frequency_df1)
#     sig_clstr2_wghtd_score = weighted_annotation(sig_dist2, frequency_df2)
#
#     """ Visualize (Approach 4) """
#     plot_figures(sig_clstr1_wghtd_score, sig_clstr2_wghtd_score,
#                  y_lable='sigmoid_weighted_annotation_score_per_cell_k_' + str(k_val))
#
#     # the weighted distance for a cluster
#     print("Approach 4: sigmoid weighted annotation :", cluster1, " ", common_cluster, sig_clstr1_wghtd_score.mean(),
#           " k: ", k_val)
#     print("Approach 4: sigmoid weighted annotation :", cluster2, " ", common_cluster, sig_clstr2_wghtd_score.mean(),
#           " k: ", k_val)
#
#
# """ Approach 5: entropy of sigmoid weighted annotation """
# entropy_wgtd_antn_clstr1 = compute_entropy_mat(sig_dist1*frequency_df1/50)
# entropy_wgtd_antn_clstr2 = compute_entropy_mat(sig_dist2*frequency_df2/50)
#
# print("Approach 5: entropy weighted annotation :", cluster1, " ", common_cluster, entropy_wgtd_antn_clstr1.mean())
# print("Approach 5: entropy weighted annotation :", cluster2, " ", common_cluster, entropy_wgtd_antn_clstr2.mean())
#
# """ Visualize (Approach 5) """
