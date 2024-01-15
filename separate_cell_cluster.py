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
# from scimilarity.utils import lognorm_counts
# from scimilarity import CellAnnotation, align_dataset
import utils_AT
import numpy as np
import pandas as pd
import csv

""""  ###### VARIABLES ########  """

annotation_path = '../models/annotation_model_v1'
data_path = '../data/KPTracer-Data/expression/adata_processed_comp_SCANTN.h5ad'

""""  ###### VARIABLES ########  """

# Model read, data read
sc.set_figure_params(dpi_save=800)
warnings.filterwarnings('ignore')
# ca = CellAnnotation(model_path=annotation_path)
adams2 = sc.read(data_path)

tumor_names = adams2.obs['Tumor']
tumor_names_arr = np.array(tumor_names.unique())
KP_tumor_list = []
for tumor in tumor_names_arr:
     try:
          if tumor.split("_")[1] != 'NT':
               KP_tumor_list.append(tumor)
     except:
          print()

# adams = adams2[adams2.obs.iloc[:, -1].isin(['Mesenchymal-1'])]
adams = adams2[adams2.obs.iloc[:, 6].isin(KP_tumor_list)]

tumor_names = adams2.obs['Tumor']
tumor_names_arr = np.array(tumor_names.unique())
mouseSpcfc_tumor_list = []
for tumor in tumor_names_arr:
     try:
          if tumor.split("_")[0] == '3457':
               mouseSpcfc_tumor_list.append(tumor)
     except:
          print()

adams = adams[adams.obs.iloc[:, 6].isin(mouseSpcfc_tumor_list)]
subset_adams1 = adams[adams.obs['Cluster-Name'].isin(['Mesenchymal-1'])]
subset_adams2 = adams[adams.obs['Cluster-Name'].isin(['Early gastric'])]
#
# # save
subset_adams1.write_h5ad('../data/yang_Mesenchymal_1_3457.h5ad')
subset_adams2.write_h5ad('../data/yang_Early gastric_3457.h5ad')


print("done")
