# Environment settings
import scanpy as sc
import warnings
from scimilarity import CellAnnotation
import utils_AT

""""  ###### VARIABLES ########  """

# Select visualization options
ALL = 0
clstr = 0
uncnstrnd_prdctn = 0
cnstrnd_prdctn = 0
min_dist = 1
gene_score = 1

annotation_path = '../models/annotation_model_v1'
# data_path = '../data/GSE136831_subsample.h5ad'
# data_path = '../data/KPTracer-Data/expression/adata_processed.combined.h5ad' #(AT)
# data_path = '../data/KPTracer-Data/expression/adata_processed_comp.h5ad' #(AT)
data_path = '../data/MSK_tumor_data/adata_combined_mnnc_010920_SCANTN2.h5ad'
# data_path = '../data/MSK_tumor_data/adata.combined.mnnc.010920.h5ad'
imag_name1 = utils_AT.create_image_name("_ MSK_tumorClstrCOMP_")
imag_name2 = utils_AT.create_image_name("_ MSK_tumorConstrndAntnCOMP_")
imag_name3 = utils_AT.create_image_name("_ MSK_tumorConstrneAntnCOMP_")
imag_name4 = utils_AT.create_image_name("_ MSK_tumorMinDistCOMP_")
imag_name5 = utils_AT.create_image_name("_ MSK_tumorGeneScoreCOMP_")
imag_name6 = utils_AT.create_image_name("_ MSK_tumorGeneScoreDscrtCOMP_")

""""  ###### VARIABLES ########  """

# Model read, data read
sc.set_figure_params(dpi_save=800)
warnings.filterwarnings('ignore')
ca = CellAnnotation(model_path=annotation_path)
adams = sc.read(data_path)


# Visualization
# if ALL or clstr:
#     sc.pl.umap(adams, color='cell_type_general', legend_fontsize=8, save=utils_AT.create_image_name('MSK_CG'), legend_loc='on data')
#     sc.pl.umap(adams, color='cell_type_fine', legend_fontsize=8, save=utils_AT.create_image_name('MSK_CF'), legend_loc='on data')
#     sc.pl.umap(adams, color='tissue', legend_fontsize=8, save=utils_AT.create_image_name('MSK_tissue'), legend_loc='on data')
#     sc.pl.umap(adams, color='histo', legend_fontsize=8, save=utils_AT.create_image_name('MSK_histo'), legend_loc='on data')
#     sc.pl.umap(adams, color='cell_type_coarse', legend_fontsize=8, save=utils_AT.create_image_name('MSK_CC'), legend_loc='on data')
#     sc.pl.umap(adams, color='cell_type_med', legend_fontsize=8, save=utils_AT.create_image_name('MSK_CM'),
#                legend_loc='on data')
#
#     sc.pl.umap(adams, color='cell_type_general', legend_fontsize=8, save=utils_AT.create_image_name('MSK_CG_legSide'))
#     sc.pl.umap(adams, color='cell_type_fine', legend_fontsize=8, save=utils_AT.create_image_name('MSK_CF_legSide'))
#     sc.pl.umap(adams, color='tissue', legend_fontsize=8, save=utils_AT.create_image_name('MSK_tissue_legSide'))
#     sc.pl.umap(adams, color='histo', legend_fontsize=8, save=utils_AT.create_image_name('MSK_histo_legSide'))
#     sc.pl.umap(adams, color='cell_type_coarse', legend_fontsize=8, save=utils_AT.create_image_name('MSK_CC_legSide'))
#     sc.pl.umap(adams, color='cell_type_med', legend_fontsize=8, save=utils_AT.create_image_name('MSK_CM_legSide'))

""" CLUSTERTING """


""" UNCONSTRAINED ANNOTATION """

# - nn_idxs: indicies of cells in the SCimilarity reference.
# - nn_dists: the minimum distance within k=50 nearest neighbors.
# - nn_stats: a dataframe containing useful metrics
#           such as (distribution of celltypes in k=50 nearest neighbors.

# Processing
predictions, nn_idxs, nn_dists, nn_stats = ca.get_predictions_kNN(adams.obsm['X_scimilarity'])
adams.obs['predictions_unconstrained'] = predictions.values

# puts in scimilarity stats as h5ad form
for ix in range(len(nn_stats.columns)):
    name = "sc_" + nn_stats.columns[ix]
    adams.obs[name] = nn_stats[nn_stats.columns[ix]].values


# (AT) get the distribution of predicted annotation
celltype_counts = adams.obs.predictions_unconstrained.value_counts()
well_represented_celltypes = celltype_counts[celltype_counts>20].index

# top 19 unconstrained predictions
celltype_counts = adams.obs.predictions_unconstrained.value_counts()
target_celltypes = celltype_counts[:19]

# adams.write('../data/MSK_tumor_data/adata_combined_mnnc_010920_SCANTN.h5ad')

# Visualization
if ALL or uncnstrnd_prdctn:
    sc.pl.umap(adams[adams.obs.predictions_unconstrained.isin(well_represented_celltypes)],
           color='predictions_unconstrained',
           legend_fontsize=8, save=imag_name2, legend_loc = "on data")
    sc.pl.umap(adams[adams.obs.predictions_unconstrained.isin(well_represented_celltypes)],
               color='predictions_unconstrained',
               legend_fontsize=8, save='legSide_'+imag_name2)


""" UNCONSTRAINED ANNOTATION """


""" CONSTRAINED ANNOTATION """

# target_celltypes = ['secretory cell', 'glandular epithelial cell', 'myofibroblast cell', 'epithelial cell', 'cardiac muscle cell', 'native cell',
#                     'luminal epithelial cell of mammary gland', 'epithelial cell of proximal tubule', 'type II pneumocyte', 'basal cell',
#                     'mesenchymal cell', 'kidney collecting duct principal cell']
# # target_celltypes = ['alveolar macrophage', 'macrophage', 'natural killer cell', 'ciliated cell', 'mature NK T cell',
#                     'B cell', 'fibroblast', 'classical monocyte', 'type II pneumocyte', 'endothelial cell of vascular tree',
#                     'club cell', 'endothelial cell of lymphatic vessel', 'CD8-positive, alpha-beta T cell',
#                     'respiratory basal cell', 'mast cell', 'type I pneumocyte', 'secretory cell', 'CD4-positive, alpha-beta T cell',
#                     'lung macrophage', 'plasma cell', 'basal cell', 'non-classical monocyte', 'plasmacytoid dendritic cell',
#                     'lung ciliated cell', 'vascular associated smooth muscle cell', 'conventional dendritic cell',
#                     'goblet cell', 'smooth muscle cell', 'pericyte', 'regulatory T cell', 'myofibroblast cell',
#                     'neuroendocrine cell', 'pulmonary ionocyte']
# Constrained classfication (reduces the noise)
target_celltypes = adams.obs['predictions_unconstrained'].value_counts()[:20].index.tolist()
ca.safelist_celltypes(target_celltypes)
adams = ca.annotate_dataset(adams, skip_preprocessing=True) # we already pre-processed the data

# Visualization
if ALL or cnstrnd_prdctn:
    sc.pl.umap(adams, color='celltype_hint', legend_fontsize=5, save=imag_name3, legend_loc = "on data")
    sc.pl.umap(adams, color='celltype_hint', legend_fontsize=5, save='legSide_'+imag_name3)

""" CONSTRAINED ANNOTATION """

""" Min Distance """

#Cell annotation computes QC metrics (min_dist), The greater min_dist the less confidence we have in the model’s prediction.
adams.obs['scimilarity_min_dist'] = nn_stats.min_dist.values

adams.write('../data/MSK_tumor_data/adata_combined_mnnc_010920_SCANTN2.h5ad')


# Visualization
if ALL or min_dist:
    sc.pl.umap(adams, color='scimilarity_min_dist', vmax=.1, save=imag_name4)

""" Min Distance """

""" Gene Score Signature """

# we input a query profile representing the cell state of interest, querycell state
#fm_basic_signature = ['SPP1', 'TREM2', 'GPNMB', 'MMP9', 'CHIT1', 'CHI3L1']
# fm_basic_signature = (utils_AT.read_csv_column
#                       ("../files/common_adata_targetGene_2023_9_20_20_54_19.csv", "Common_gene_name"))
# sc.tl.score_genes(adams, fm_basic_signature)
# if ALL or gene_score:
#     sc.pl.umap(adams, color='score', save=imag_name5)
#
# adams.write('../data/KPTracer-Data/expression/adata_processed_comp_SCANTN2.h5ad')
# print()
# # make a subset of anndata which gene score more than 0.5
# threshold = 0.5
# values = adams.obs['score'].values
# colors = np.where(values > threshold, 'blue', 'gray')
# adams.obs['gene_score_dscrt'] = colors
#
# # we are creating a subset just with score more than threlhold
# subset_ftr = ['blue']
# adamsUPgnscr = adams[adams.obs.iloc[:, -1].isin(subset_ftr)]
#
#
# yang_clstrs = adamsUPgnscr.obs['Cluster-Name'].unique()
# for clstr in yang_clstrs:
#     numOfCells = len(adamsUPgnscr[adamsUPgnscr.obs['Cluster-Name'] == clstr])
#     ratio = round(100 * (numOfCells / adamsUPgnscr.shape[0]), 2)
#     print(clstr, " ", ratio)

""" Gene Score Signature """

print("done")
