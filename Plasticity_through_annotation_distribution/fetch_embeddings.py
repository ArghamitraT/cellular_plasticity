""" this files fetches the unique cell types from train_meta.csv file and fetches their corresponding embeddings """


import os
import json
import numpy as np

# Define the data paths
datapath = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/index_values_by_type.json')
embedding_path = os.path.join(os.getcwd(), '../../models/query_model_v1/train_embedding.npy')

# Load the index values and the embeddings
with open(datapath, 'r') as fp:
    index_values_by_type = json.load(fp)

embeddings = np.load(embedding_path)

# Create a directory to save cell type embeddings
output_dir = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/cell_type_embeddings')
os.makedirs(output_dir, exist_ok=True)

# Loop through each cell type and its indices
for cell_type, indices in index_values_by_type.items():
    print(f"Processing cell type: {cell_type}")

    # Convert indices to integers
    cell_indices = [int(index) for index in indices]

    # Fetch the rows from the embeddings with the cell indices
    cell_embeddings = embeddings[cell_indices]

    # Define the filename to save the embeddings
    cell_type = cell_type.replace(" ", "_")
    output_filename = f"{cell_type}_embeddings.npy"
    output_filepath = os.path.join(output_dir, output_filename)

    # Save the embeddings as a .npy file
    np.save(output_filepath, cell_embeddings)

    print(f"Saved embeddings for {cell_type} to {output_filepath}")
