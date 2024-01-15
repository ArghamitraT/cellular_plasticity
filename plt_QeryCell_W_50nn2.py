"""
This file plots selected query cells along its 50 nearest neighbors.
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
from collections import Counter


data_path = '../data/MSK_tumor_data/adata_combined_mnnc_010920_SCANTN2.h5ad'
adata = sc.read(data_path)

""" choosing a part of UMAP """
sc.pl.umap(adata, color='scimilarity_min_dist', vmax=.1, show=False)
#center_x, center_y = 12.562, -3.798  # low confidence
center_x, center_y = -4.871, 1.543     # High confidence
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
    min_dist_data = adata.obs['sc_all_nn_indices'][idx]
    nn_names = adata.obs['sc_all_nn_names'][idx]  # Get nn names for each query cell
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
        # nn_indices.extend([int(i) for i in decoded_data])       # Convert each element in the list to an integer
        # nn_names_list.extend([i for i in decoded_names])
        nn_indices.extend([decoded_data[0]])  # Convert each element in the list to an integer
        nn_names_list.extend([decoded_names[0]])
    else:
        nn_indices.append(int(decoded_data))        # Convert a single element to an integer and append


""" Fetch 128D embeddings """
all_embeddings = np.load('../models/query_model_v1/train_embedding.npy')
nn_embeddings = all_embeddings[nn_indices, :]
query_embeddings = adata.obsm['X_scimilarity'][query_cell_indices]


""" Get all UMAPS """
# Combine the query embeddings with the nn embeddings
combined_embeddings = np.vstack([query_embeddings, nn_embeddings])
# Perform UMAP reduction on the combined subset
reducer = umap.UMAP()
combined_umap = reducer.fit_transform(combined_embeddings)


""" Recalculate the UMAP only with query cells and their 50 nearest neighbors and then plot. """
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
# Plotting the new UMAP coordinates with labels
plt.figure(figsize=(14, 10))  # Adjust the figure size as needed
scatter = sns.scatterplot(
    x=combined_umap[:, 0],
    y=combined_umap[:, 1],
    hue=labels,
    palette=palette,
    legend='brief',  # Use 'brief' or 'full' as per your preference
    s=10  # Size of the points
)
# Place the legend to the side
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
# Set title and labels
plt.title('UMAP visualization with Query Cells and NN Embeddings')
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
# Adjust the plot area to make room for the legend
plt.tight_layout(rect=[0, 0, 0.75, 1])
# Save the figure with a suitable name
plt.savefig('figures/' + utils_AT.create_image_name('recalculatedUMAPWqcWnn'))
# Show the plot
plt.show()



