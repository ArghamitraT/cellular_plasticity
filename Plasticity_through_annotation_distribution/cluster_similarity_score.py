import os
import numpy as np
import pandas as pd
import scanpy as sc
import utils_AT as utils
from scipy.stats import multivariate_normal

data_path = '../../data/KPTracer-Data/expression/adata_processed_comp_SCANTN2.h5ad'
adams_comp = sc.read(data_path)

# Filter the data
adams1 = adams_comp[adams_comp.obs['Cluster-Name'].isin(['Endoderm-like'])]
adams1 = adams1[adams1.obs['celltype_hint'].isin(['secretory cell'])]

# Suppose `query_cell` is the new data point for which you want to calculate the log-likelihood
query_cell = adams1.obsm['X_scimilarity']

# Prepare the frequency matrix
freq_mat = utils.make_frequency_matrix(adams1.obs['sc_hits'])
freq_mat.columns = freq_mat.columns.str.replace(' ', '_')
cell_types_to_consider = freq_mat.columns.tolist()