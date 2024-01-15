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
# from bokeh.plotting import figure, show, output_file
# from bokeh.models import HoverTool, ColumnDataSource
# from bokeh.transform import linear_cmap
# from bokeh.palettes import Viridis256


data_path = '../data/MSK_tumor_data/adata_combined_mnnc_010920_SCANTN2.h5ad'
adata = sc.read(data_path)

""" Interactive UMAP """
# Assuming 'adata' is your AnnData object and UMAP has been computed on it
# umap_coords = adata.obsm['X_umap']
#
#
# # Assuming 'umap_coords' are already extracted from 'adata'
# source = ColumnDataSource(data=dict(
#     x=umap_coords[:, 0],
#     y=umap_coords[:, 1],
#     desc=[f'Cell {i}' for i in range(umap_coords.shape[0])]
# ))
#
# output_file("umap_basic_plot.html")
#
# tools = "pan,wheel_zoom,box_select,lasso_select,reset"
# p = figure(tools=tools, title="UMAP Plot", toolbar_location="above")
#
# p.circle('x', 'y', source=source, size=5)
#
# hover = HoverTool()
# hover.tooltips = [("Index", "@desc"), ("(x,y)", "($x, $y)")]
# p.add_tools(hover)
#
# show(p)

""" choosing a part of UMAP """
sc.pl.umap(adata, color='scimilarity_min_dist', vmax=.1, show=False)
#center_x, center_y = 12.562, -3.798  # Adjust according to your plot
center_x, center_y = -4.871, 1.543
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
nn_names = []
for idx in query_cell_indices:
    min_dist_data = adata.obs['sc_min_dist_index'][idx]
    if isinstance(min_dist_data, str):       # Check if the data is a string (and hence, possibly a JSON string)
        try:
            decoded_data = json.loads(min_dist_data)        # Try to decode the JSON string
        except json.JSONDecodeError:
            print(f"Failed to decode JSON for index {idx}: {min_dist_data}")
            continue
    else:
        decoded_data = min_dist_data        # If it's not a string, use the data directly
    if isinstance(decoded_data, list):      # Now process the decoded data
        nn_indices.extend([int(i) for i in decoded_data])       # Convert each element in the list to an integer
    else:
        nn_indices.append(int(decoded_data))        # Convert a single element to an integer and append


""" Fetch 128D embeddings """
all_embeddings = np.load('../models/query_model_v1/train_embedding.npy')
nn_embeddings = all_embeddings[nn_indices, :]


""" Get all UMAPS """
reducer = umap.UMAP()
nn_umap_coords = reducer.fit_transform(nn_embeddings)       # If 'nn_embeddings' are the embeddings of nearest neighbors, first transform them using UMAP and then plot
query_umap_coords = adata.obsm['X_umap'][query_cell_indices, :]     # Highlight query cells
umap_coords = adata.obsm['X_umap']      # Get the entire UMAP


""" Plotting UMAP with all cells """
fig, ax = plt.subplots(figsize=(12, 8))
sns.scatterplot(ax=ax, x=umap_coords[:, 0], y=umap_coords[:, 1], color="grey", label="All Cells")
sns.scatterplot(ax=ax, x=nn_umap_coords[:, 0], y=nn_umap_coords[:, 1],
                color="blue", label="Nearest Neighbors")        # Highlight nearest neighbors
sns.scatterplot(ax=ax, x=query_umap_coords[:, 0], y=query_umap_coords[:, 1], color="red", label="Query Cells")
ax.set_title("UMAP Plot with Query Cells and their closest neighbor")
ax.set_xlim(-20, 40)
ax.set_ylim(-20, 40)
plt.legend()
plt.savefig('figures/'+utils_AT.create_image_name('QCvsRefCellsWallCells'))
plt.show()


""" Plotting only with query cells and its nn """
fig, ax = plt.subplots(figsize=(12, 8))
sns.scatterplot(ax=ax, x=nn_umap_coords[:, 0], y=nn_umap_coords[:, 1],
                color="blue", label="Nearest Neighbors")        # Highlight nearest neighbors
sns.scatterplot(ax=ax, x=query_umap_coords[:, 0], y=query_umap_coords[:, 1], color="red", label="Query Cells")
plt.title('UMAP visualization of Query Cells and NN Embeddings')
ax.set_xlim(-20, 40)
ax.set_ylim(-20, 40)
plt.legend()
plt.savefig('figures/'+utils_AT.create_image_name('QCvsRefCells'))
plt.show()


""" Recalculate the UMAP with all the cells and then plot. """
# Load the 128-dimensional embeddings from your anndata object
existing_embeddings = adata.obsm['X_scimilarity']
# Combine the existing embeddings with the additional cells' data
combined_embeddings = np.vstack([existing_embeddings, nn_embeddings])
# Perform UMAP reduction on the combined dataset
reducer = umap.UMAP()
combined_umap = reducer.fit_transform(combined_embeddings)
# Create labels for the plot
# 0 for all original cells, 1 for additional cells, 2 for query cells
labels = np.zeros(combined_embeddings.shape[0])
labels[-nn_embeddings.shape[0]:] = 1  # Label for additional cells
labels[query_cell_indices] = 2  # Label for query cells
# Plotting the new UMAP coordinates with labels
plt.figure(figsize=(12, 8))
sns.scatterplot(x=combined_umap[:, 0], y=combined_umap[:, 1], hue=labels, palette=['grey', 'blue', 'red'], legend='full', s=10)
plt.title('Combined UMAP visualization with Original, Additional, and Query Cells')
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.savefig('figures/'+utils_AT.create_image_name('recalculatedUMAP'))
plt.show()

""" Recalculate the UMAP only with query cells and their nearest neighbors and then plot. """
# Extract the embeddings for the query cells
query_embeddings = existing_embeddings[query_cell_indices]
# Combine the query embeddings with the nn embeddings
combined_embeddings = np.vstack([query_embeddings, nn_embeddings])
# Perform UMAP reduction on the combined subset
reducer = umap.UMAP()
combined_umap = reducer.fit_transform(combined_embeddings)
# Create labels for the plot: 0 for query cells, 1 for nn embeddings
labels = np.zeros(combined_embeddings.shape[0])
labels[len(query_embeddings):] = 1  # Label for nn embeddings
# Plotting the new UMAP coordinates with labels
plt.figure(figsize=(12, 8))
sns.scatterplot(x=combined_umap[:, 0], y=combined_umap[:, 1], hue=labels, palette=['red', 'blue'], legend='full', s=10)
plt.title('UMAP visualization with Query Cells and NN Embeddings')
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.savefig('figures/'+utils_AT.create_image_name('recalculatedUMAPWqcWnn'))
plt.show()