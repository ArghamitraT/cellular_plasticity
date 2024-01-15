"""
This file fits a UMAP on given adata, then projects new points
"""

import pandas as pd
import utils_AT
import json
import numpy as np
import seaborn as sns
import umap
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pickle


data_path = '../data/MSK_tumor_data/adata_combined_mnnc_010920_SCANTN2.h5ad'
adata = sc.read(data_path)

""" choosing a part of UMAP """
sc.pl.umap(adata, color='scimilarity_min_dist', vmax=.1, show=False)
center_x, center_y = 12.562, -3.798  # low confidence
# center_x, center_y = -4.871, 1.543     # High confidence
umap_coords = adata.obsm['X_umap']
distances = np.sqrt((umap_coords[:, 0] - center_x)**2 + (umap_coords[:, 1] - center_y)**2)
num_cells = 500
selected_cells_indices = np.argsort(distances)[:num_cells]
radius = np.max(distances[selected_cells_indices])
circle = patches.Circle((center_x, center_y), radius, edgecolor='red', facecolor='none', linewidth=2)
ax = plt.gca()
ax.add_patch(circle)
plt.title('UMAP Plot with Highlighted Region')
plt.savefig('figures/'+utils_AT.create_image_name('selected_region'))
plt.show()


""" Get Query cell and their closest neighbor indices """
query_cell_indices = selected_cells_indices.tolist()
nn_indices = []
nn_names_list = []
for idx in query_cell_indices:
    # min_dist_data = adata.obs['sc_all_nn_indices'][idx]
    # nn_names = adata.obs['sc_all_nn_names'][idx]  # Get nn names for each query cell
    min_dist_data = adata.obs['sc_all_nn_indices'].iloc[idx]
    nn_names = adata.obs['sc_all_nn_names'].iloc[idx]  # Get nn names for each query cell
    if isinstance(min_dist_data, str):       # Check if the data is a string (and hence, possibly a JSON string)
        try:
            decoded_data = json.loads(min_dist_data)        # Try to decode the JSON string
            decoded_names = json.loads(nn_names)
        except json.JSONDecodeError:
            print(f"Failed to decode JSON for index {idx}: {min_dist_data}")
            continue
    else:
        decoded_data = min_dist_data        # If it's not a string, use the data directly
    if isinstance(decoded_data, list):      # Now process the decoded data
        nn_indices.extend([int(i) for i in decoded_data])       # Convert each element in the list to an integer
        nn_names_list.extend([i for i in decoded_names])
        # nn_indices.extend([decoded_data[0]])  # Convert each element in the list to an integer
        # nn_names_list.extend([decoded_names[0]])
    else:
        nn_indices.append(int(decoded_data))        # Convert a single element to an integer and append


""" Fetch 128D embeddings """
all_embeddings = np.load('../models/query_model_v1/train_embedding.npy')
nn_embeddings = all_embeddings[nn_indices, :]
all_query_embeddings = adata.obsm['X_scimilarity']


""" fit the training and test """
# fit all QC and save
# trans = umap.UMAP(random_state=42).fit(all_query_embeddings)
# with open('umap_model.pkl', 'wb') as output_file:
#     pickle.dump(trans, output_file)
# Load the model from the file; Transforming and plotting the test set embedding
with open('umap_model.pkl', 'rb') as input_file:
    trans = pickle.load(input_file)
test_embedding = trans.transform(nn_embeddings)

# """ Plot the all QC and reference nn individually"""
# # Plotting the training set embedding
# plt.figure(figsize=(10, 8))
# # Plot all cells in grey
# sns.scatterplot(x=trans.embedding_[:, 0], y=trans.embedding_[:, 1], color="grey", alpha=0.5)
# # Overlay query cells in red
# sns.scatterplot(x=trans.embedding_[query_cell_indices, 0], y=trans.embedding_[query_cell_indices, 1], color="red", label='Query Cells')
# plt.title('Embedding of the training set by UMAP', fontsize=24)
# plt.legend()
# plt.savefig('figures/' + utils_AT.create_image_name('AllQCell'))
# plt.show()
# # plot the reference nn
# plt.figure(figsize=(10, 8))
# sns.scatterplot(x=test_embedding[:, 0], y=test_embedding[:, 1], s=5, color="blue", label='NN Embeddings')
# plt.title('Embedding of the test set by UMAP', fontsize=24)
# # plt.title('Embedding of the test set (50 nn) by UMAP', fontsize=24)
# plt.legend()
# plt.savefig('figures/' + utils_AT.create_image_name('projectedRefCell'))
# plt.show()
#
#
# """ plot the training and test UMAP together """
# # Create a color list for the combined embeddings
# combined_colors = ['grey'] * len(trans.embedding_)
# # Update colors for query cells to red
# for idx in query_cell_indices:
#     combined_colors[idx] = 'red'
# # Extend the color list with blue for test embeddings
# combined_colors.extend(['blue'] * len(test_embedding))
# # Combine the trans embedding with the test embedding
# combined_embedding = np.vstack([trans.embedding_, test_embedding])
# # Plotting the combined UMAP embeddings
# plt.figure(figsize=(10, 8))
# # Create a DataFrame for the combined embeddings and colors
# embedding_df = pd.DataFrame(combined_embedding, columns=['UMAP_1', 'UMAP_2'])
# embedding_df['Color'] = combined_colors
# # Plot using seaborn
# sns.scatterplot(data=embedding_df, x='UMAP_1', y='UMAP_2', hue='Color', palette=['grey', 'red', 'blue'], alpha=0.5)
# plt.title('Combined UMAP Embedding', fontsize=24)
# plt.legend()
# plt.savefig('figures/' + utils_AT.create_image_name('CombinedUMAP'))
# plt.show()


""" just plot the selected query cells and nn """
# Combine the relevant part of trans embedding (query cells) with the test embedding
query_embedding_subset = trans.embedding_[query_cell_indices]
combined_embedding = np.vstack([query_embedding_subset, test_embedding])
# # Create a color list for the combined embeddings. Colors for query cells (red) and test embeddings (blue)
# combined_colors = ['red'] * len(query_embedding_subset) + ['blue'] * len(test_embedding)
# # Plotting the combined UMAP embeddings
# plt.figure(figsize=(10, 8))
# # Create a DataFrame for the combined embeddings and colors
# embedding_df = pd.DataFrame(combined_embedding, columns=['UMAP_1', 'UMAP_2'])
# embedding_df['Color'] = combined_colors
# # Plot using seaborn
# sns.scatterplot(data=embedding_df, x='UMAP_1', y='UMAP_2', hue='Color', palette=['red', 'blue'], alpha=0.5)
# # plt.title('Combined UMAP Embedding (only QC and 50 nn)', fontsize=24)
# plt.title('Combined UMAP Embedding (only QC and nn)', fontsize=24)
# plt.legend()
# plt.savefig('figures/' + utils_AT.create_image_name('CombinedUMAP'))
# plt.show()


""" Plot reference cells with cell types """
# Prepare labels for coloring
labels = ['Query'] * len(query_cell_indices)  # Label query cells as "Query"
# For each nearest neighbor cell, use their cell type for labeling
labels.extend(nn_names_list)  # Append these labels to the label list
# Create a set of unique labels excluding 'Query'
unique_labels = set(labels) - {'Query'}
# Define a color palette for non-query cells
# You can choose a specific palette or use seaborn's color palette functions
non_query_palette = sns.color_palette("husl", len(unique_labels))
# Create a complete palette with red for 'Query' and other colors for non-query cells
palette = {'Query': 'red'}
palette.update(dict(zip(unique_labels, non_query_palette)))
plt.figure(figsize=(14, 10))  # Adjust the figure size as needed
scatter = sns.scatterplot(
    x=combined_embedding[:, 0],
    y=combined_embedding[:, 1],
    hue=labels,
    palette=palette,
    legend='brief',  # Use 'brief' or 'full' as per your preference
    s=10  # Size of the points
)
plt.savefig('figures/' + utils_AT.create_image_name('CombinedUMAPwCellType'))
plt.show()