# Environment settings
import scanpy as sc
import warnings
from scimilarity.utils import lognorm_counts
from scimilarity import CellAnnotation, align_dataset

###### VARIABLES ########
annotation_path = '../../models/annotation_model_v1'
#data_path = '../data/GSE136831_subsample.h5ad'
data_path = '../../data/KPTracer-Data/expression/adata_processed.combined.h5ad'  #(AT)
import utils_AT         #(AT) image names
imag_name1 = utils_AT.create_image_name("_DyangDataClstr_")
imag_name2 = utils_AT.create_image_name("_DyangDataAntn_")
imag_name3 = utils_AT.create_image_name("_DyangDataConstrneClss_")
imag_name4 = utils_AT.create_image_name("_DyangDataMinDist_")
imag_name5 = utils_AT.create_image_name("_DyangDataQueryClStat_")
###### VARIABLES ########

sc.set_figure_params(dpi=100)
warnings.filterwarnings('ignore')
ca = CellAnnotation(model_path=annotation_path)
adams = sc.read(data_path)

#(AT) to match the target gene, all name of the gene needs to be uppercase
for indx_num in range(len(adams.var.index)):
     adams.var.index.values[indx_num] = adams.var.index[indx_num].upper()

#data processing
adams = align_dataset(adams, ca.gene_order)
adams = lognorm_counts(adams)
adams.obsm['X_scimilarity'] = ca.get_embeddings(adams.X)

#clustering and visualization
sc.pp.neighbors(adams, use_rep='X_scimilarity')
sc.tl.umap(adams)
#sc.pl.umap(adams, color='celltype_raw', legend_fontsize=8, save=imag_name1)    #(uncomment LATER (AT))

#Unconstrained annotation
"""
- nn_idxs: indicies of cells in the SCimilarity reference. 
- nn_dists: the minimum distance within k=50 nearest neighbors. 
- nn_stats: a dataframe containing useful metrics 
          such as (distribution of celltypes in k=50 nearest neighbors.
"""
# (AT) in the nn_dists we have min_dist
predictions, nn_idxs, nn_dists, nn_stats = ca.get_predictions_kNN(adams.obsm['X_scimilarity'])
adams.obs['predictions_unconstrained'] = predictions.values
celltype_counts = adams.obs.predictions_unconstrained.value_counts()
well_represented_celltypes = celltype_counts[celltype_counts>20].index

#(AT)trying to find out Dyang's cell distribution
#Dyang_cell_dist = adams.obs.iloc[:, -1].value_counts()
#utils_AT.bar_plt_distribution(Dyang_cell_dist, title='DYang cell Distribution', image_name='_DYangCellDist_')

#(AT) to get a distribution on # of cell types
#max_cellcount = celltype_counts[:15]
#utils_AT.bar_plt_distribution(max_cellcount)

# Annotation visualization (uncomment LATER (AT))
# sc.pl.umap(adams[adams.obs.predictions_unconstrained.isin(well_represented_celltypes)],
#            color='predictions_unconstrained',
#            legend_fontsize=8, save=imag_name2)

#only show the cell annotation for high plasticity cells
subset_ftr = ['High plasticity']
subset_adams = adams[adams.obs.iloc[:, -2].isin(subset_ftr)]

sc.pl.umap(subset_adams[subset_adams.obs.predictions_unconstrained.isin(well_represented_celltypes)],
           color='predictions_unconstrained',
           legend_fontsize=5, save=imag_name2, legend_loc = "on data")

#Constrained classfication (reduces the noise) (uncomment LATER (AT))
target_celltypes = ['secretory cell', 'glandular epithelial cell', 'myofibroblast cell', 'epithelial cell', 'cardiac muscle cell', 'native cell',
                    'luminal epithelial cell of mammary gland', 'epithelial cell of proximal tubule', 'type II pneumocyte', 'basal cell']

# target_celltypes = ['alveolar macrophage', 'macrophage', 'natural killer cell', 'ciliated cell', 'mature NK T cell',
#                     'B cell', 'fibroblast', 'classical monocyte', 'type II pneumocyte', 'endothelial cell of vascular tree',
#                     'club cell', 'endothelial cell of lymphatic vessel', 'CD8-positive, alpha-beta T cell',
#                     'respiratory basal cell', 'mast cell', 'type I pneumocyte', 'secretory cell', 'CD4-positive, alpha-beta T cell',
#                     'lung macrophage', 'plasma cell', 'basal cell', 'non-classical monocyte', 'plasmacytoid dendritic cell',
#                     'lung ciliated cell', 'vascular associated smooth muscle cell', 'conventional dendritic cell',
#                     'goblet cell', 'smooth muscle cell', 'pericyte', 'regulatory T cell', 'myofibroblast cell',
#                     'neuroendocrine cell', 'pulmonary ionocyte']
ca.safelist_celltypes(target_celltypes)
adams = ca.annotate_dataset(adams, skip_preprocessing=True) # we already pre-processed the data

sc.pl.umap(adams, color='celltype_hint', legend_fontsize=5, save=imag_name3, legend_loc = "on data")

#Cell annotation computes QC metrics (min_dist)
#The greater min_dist the less confidence we have in the model’s prediction.
#DID NOT WORK; look into it later (AT)
#sc.pl.umap(adams, color='min_dist', vmax=.1, save=imag_name4)

#we input a query profile representing the cell state of interest, querycell state
fm_basic_signature = ['SPP1', 'TREM2', 'GPNMB', 'MMP9', 'CHIT1', 'CHI3L1']
sc.tl.score_genes(adams, fm_basic_signature)
sc.pl.umap(adams, color='score', save=imag_name5)

print("done")
