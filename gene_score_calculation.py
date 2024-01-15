""" For MSK lung tumor data, this files finds the gene score
(for aligned genes with scimilarity) of mentioned cell types """

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import scipy.sparse

# Load the single-cell data (assuming 'adata' is your Anndata object)
data_path = '../data/MSK_tumor_data/adata_combined_mnnc_010920_SCANTN2.h5ad'
adata = sc.read(data_path)  # Load your data here
type1 = 'SCLC-A'
type2 = 'T cell'

# Load the CSV file containing important genes
important_genes = pd.read_csv('../files/MSK_lung_tumor_aligned_genes.csv')
target_genes = important_genes['Gene'].tolist()  # Assuming gene names are in a column named 'Gene'

# Filter the data for the two cell types
cell_type_1 = adata[adata.obs['cell_type_med'] == type1]
cell_type_2 = adata[adata.obs['cell_type_med'] == type2]

# Function to calculate average gene score
# def average_gene_score(cell_data, genes):
#     return cell_data[:, genes].X.mean(axis=1)

def average_gene_score(cell_data, genes):
    # Extract the gene expression data
    gene_expression = cell_data[:, genes].X

    # If the data is in a sparse matrix format, convert it to a dense array
    if scipy.sparse.issparse(gene_expression):
        gene_expression = gene_expression.toarray()

    # Calculate the mean across genes for each cell
    avg_scores = np.mean(gene_expression, axis=1)

    # Convert avg_scores to a 1-dimensional array if it's not already
    return np.ravel(avg_scores)


# Calculate average gene score for each cell type
avg_score_type_1 = average_gene_score(cell_type_1, target_genes)
avg_score_type_2 = average_gene_score(cell_type_2, target_genes)

# Prepare DataFrame for plotting
df_type_1 = pd.DataFrame({
    'Average Gene Score': avg_score_type_1,
    'Cell Type': [type1] * len(avg_score_type_1)
})

df_type_2 = pd.DataFrame({
    'Average Gene Score': avg_score_type_2,
    'Cell Type': [type2] * len(avg_score_type_2)
})

# Concatenate the dataframes
plot_data = pd.concat([df_type_1, df_type_2])

# Plotting the violin plot
sns.violinplot(x='Cell Type', y='Average Gene Score', data=plot_data)
plt.title('Average Gene Score by Cell Type')
plt.savefig('figures/avg_gene_score.png')
plt.show()

print()
