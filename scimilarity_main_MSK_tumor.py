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
data_path = '../data/MSK_tumor_data/adata_combined_mnnc_010920_SCANTN.h5ad'


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

adams.write('../data/MSK_tumor_data/try.h5ad')

# Visualization
if ALL or min_dist:
    sc.pl.umap(adams, color='scimilarity_min_dist', vmax=.1, save=imag_name4)

""" Min Distance """

print("done")
