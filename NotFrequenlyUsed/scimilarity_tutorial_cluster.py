# Environment settings
import scanpy as sc

#(AT)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import utils_AT

sc.set_figure_params(dpi_save=800)
import warnings
warnings.filterwarnings('ignore')
from scimilarity.utils import lognorm_counts
from scimilarity import CellAnnotation, align_dataset


# Instantiate the CellAnnotation object.
# Replace model_path with your local file path.
annotation_path = '../../models/annotation_model_v1'
ca = CellAnnotation(model_path=annotation_path)

#(AT) comment later
# data_path = '../data/GSE136831_subsample.h5ad'
# adams_atlas = sc.read(data_path)
# data_path = '../data/KPTracer-Data/expression/adata_processed.combined.h5ad' #(AT)
# adams_dyang = sc.read(data_path)
# print("done")

# Replace data_path with your local file path.
#data_path = '../data/GSE136831_subsample.h5ad'
#data_path = '../data/KPTracer-Data/expression/adata_processed.nt.h5ad' #(AT)
data_path = '../../data/KPTracer-Data/expression/adata_processed.combined.h5ad'  #(AT)
adams = sc.read(data_path)

#(AT) to match the target gene, all name of the gene needs to be uppercase
for indx_num in range(len(adams.var.index)):
     adams.var.index.values[indx_num] = adams.var.index[indx_num].upper()
sigscores = pd.read_csv(
    f"../../models/KPTracer-release-main/reproducibility/Figure6_S6/data/fitness_signature_scores.tsv",
    sep='\t', index_col = 0)
kii = np.intersect1d(sigscores.index, adams.obs_names)
adams.obs['FitnessSignature'] = np.nan
adams.obs.loc[kii, 'FitnessSignature'] = sigscores.loc[kii, 'FitnessSignature_NT']


adams = align_dataset(adams, ca.gene_order)
adams = lognorm_counts(adams)
adams.obsm['X_scimilarity'] = ca.get_embeddings(adams.X)

#(AT)
image_name = utils_AT.create_image_name("_Dyang_Atlas_")

sc.pp.neighbors(adams, use_rep='X_scimilarity')
sc.tl.umap(adams)

#sc.pl.umap(adams, color='celltype_raw', legend_fontsize=5)
sc.pl.umap(adams, color='Cluster-Name', legend_fontsize=5, save=image_name, legend_loc = "on data") #(AT)
print("done")