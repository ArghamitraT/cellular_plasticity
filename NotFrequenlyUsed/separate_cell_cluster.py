"""\ (AT)
    Given a h5ad dataset, it annotates the cells using Scimilarity and as per clusters
    it divides and save in two separate .h5ad files

    Parameters
    ----------
    adata
        Annotated data matrix.

    Returns
    -------
    Nothing.

    Saves the files in the given directory as two separate files

    Notes
    -----
    """

# Environment settings
import scanpy as sc
import warnings
from scimilarity.utils import lognorm_counts
from scimilarity import CellAnnotation, align_dataset
import utils_AT
import numpy as np
import pandas as pd
import csv

""""  ###### VARIABLES ########  """

annotation_path = '../../models/annotation_model_v1'
data_path = '../../data/KPTracer-Data/expression/adata_processed.combined.h5ad'

""""  ###### VARIABLES ########  """

# Model read, data read
sc.set_figure_params(dpi_save=800)
warnings.filterwarnings('ignore')
ca = CellAnnotation(model_path=annotation_path)
adams = sc.read(data_path)

#(AT) to match the target gene, all name of the gene needs to be uppercase
for indx_num in range(len(adams.var.index)):
     adams.var.index.values[indx_num] = adams.var.index[indx_num].upper()

#data processing
#adams_2 = align_dataset(adams, ca.gene_order)
adams = align_dataset(adams, ca.gene_order)
adams = lognorm_counts(adams)
adams.obsm['X_scimilarity'] = ca.get_embeddings(adams.X)


""" CLUSTERTING """

# Clustering
sc.pp.neighbors(adams, use_rep='X_scimilarity')
sc.tl.umap(adams)

""" CLUSTERTING """


""" UNCONSTRAINED ANNOTATION """
# - nn_idxs: indicies of cells in the SCimilarity reference.
# - nn_dists: the minimum distance within k=50 nearest neighbors.
# - nn_stats: a dataframe containing useful metrics
#           such as (distribution of celltypes in k=50 nearest neighbors.

# Processing
predictions, nn_idxs, nn_dists, nn_stats = ca.get_predictions_kNN(adams.obsm['X_scimilarity'])
adams.obs['predictions_unconstrained'] = predictions.values

""" UNCONSTRAINED ANNOTATION """

""" cluster separation for particular cluster (high plastic cells) """
# subset_ftr = ['High plasticity']
# subset_adams_HP = adams[adams.obs.iloc[:, -2].isin(subset_ftr)]
# subset_adams_HP_scrtrCell = subset_adams_HP[subset_adams_HP.obs.iloc[:, -1].isin(['secretory cell'])]
# subset_adams_HP_glndEpthCell = subset_adams_HP[subset_adams_HP.obs.iloc[:, -1].isin(['glandular epithelial cell'])]
#
# # save
# subset_adams_HP_scrtrCell.write_h5ad('../data/yang_HP_scrtrCell.h5ad')
# subset_adams_HP_glndEpthCell.write_h5ad('../data/yang_HP_glndEpthCell.h5ad')

""" cluster separation for particular cluster (Pre-EMT) """
# subset_ftr = ['Pre-EMT']
# subset_adams = adams[adams.obs.iloc[:, -2].isin(subset_ftr)]
# subset_adams1 = subset_adams[subset_adams.obs.iloc[:, -1].isin(['secretory cell'])]
# subset_adams2 = subset_adams[subset_adams.obs.iloc[:, -1].isin(['glandular epithelial cell'])]
# #
# # # save
# subset_adams1.write_h5ad('../data/yang_PrEMT_scrtrCell.h5ad')
# subset_adams2.write_h5ad('../data/yang_PrEMT_glndEpthCell.h5ad')

""" cluster separation for particular cluster (Mesen-2) """
subset_ftr = ['Mesenchymal-2']
subset_adams = adams[adams.obs.iloc[:, -2].isin(subset_ftr)]
subset_adams1 = subset_adams[subset_adams.obs.iloc[:, -1].isin(['cardiac muscle cell'])]
subset_adams2 = subset_adams[subset_adams.obs.iloc[:, -1].isin(['epithelial cell'])]

#save
subset_adams1.write_h5ad('../data/yang_Mesen2_cardiacCell.h5ad')
subset_adams2.write_h5ad('../data/yang_Mesen2_epthCell.h5ad')

""" cluster separation for particular cluster (Mesen-2 (met)) """
subset_ftr = ['Mesenchymal-2 (Met)']
subset_adams = adams[adams.obs.iloc[:, -2].isin(subset_ftr)]
subset_adams1 = subset_adams[subset_adams.obs.iloc[:, -1].isin(['cardiac muscle cell'])]
subset_adams2 = subset_adams[subset_adams.obs.iloc[:, -1].isin(['epithelial cell'])]
#
# # save
subset_adams1.write_h5ad('../data/yang_Mesen2Met_cardianCell.h5ad')
subset_adams2.write_h5ad('../data/yang_Mesen2Met_epthCell.h5ad')

""" cluster separation for particular cluster (gastric) """
subset_ftr = ['Gastric-like']
subset_adams = adams[adams.obs.iloc[:, -2].isin(subset_ftr)]
subset_adams1 = subset_adams[subset_adams.obs.iloc[:, -1].isin(['secretory cell'])]

# # save
subset_adams1.write_h5ad('../data/yang_GASTRIC_secretoryCell.h5ad')

print("done")
