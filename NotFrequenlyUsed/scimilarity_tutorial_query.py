# Environment settings
import scanpy as sc
import warnings
from scimilarity.utils import lognorm_counts
from scimilarity import CellAnnotation, align_dataset

###### VARIABLES ########
annotation_path = '../../models/annotation_model_v1'
data_path = '../../data/GSE136831_subsample.h5ad'
#data_path = '../data/KPTracer-Data/expression/adata_processed.combined.h5ad' #(AT)
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

# Annotation visualization (uncomment LATER (AT))
# sc.pl.umap(adams[adams.obs.predictions_unconstrained.isin(well_represented_celltypes)],
#            color='predictions_unconstrained',
#            legend_fontsize=8, save=imag_name2)
# sc.pl.umap(adams[adams.obs.predictions_unconstrained.isin(well_represented_celltypes)],
#            color='predictions_unconstrained',
#            legend_fontsize=2, save=imag_name2, legend_loc = "on data")

#Constrained classfication (reduces the noise) (uncomment LATER (AT))
# target_celltypes = ['alveolar macrophage', 'macrophage', 'natural killer cell', 'ciliated cell', 'mature NK T cell',
#                     'B cell', 'fibroblast', 'classical monocyte', 'type II pneumocyte', 'endothelial cell of vascular tree',
#                     'club cell', 'endothelial cell of lymphatic vessel', 'CD8-positive, alpha-beta T cell',
#                     'respiratory basal cell', 'mast cell', 'type I pneumocyte', 'secretory cell', 'CD4-positive, alpha-beta T cell',
#                     'lung macrophage', 'plasma cell', 'basal cell', 'non-classical monocyte', 'plasmacytoid dendritic cell',
#                     'lung ciliated cell', 'vascular associated smooth muscle cell', 'conventional dendritic cell',
#                     'goblet cell', 'smooth muscle cell', 'pericyte', 'regulatory T cell', 'myofibroblast cell',
#                     'neuroendocrine cell', 'pulmonary ionocyte']
# ca.safelist_celltypes(target_celltypes)
# adams = ca.annotate_dataset(adams, skip_preprocessing=True) # we already pre-processed the data
# sc.pl.umap(adams, color='celltype_hint', legend_fontsize=3, save=imag_name3)

#Cell annotation computes QC metrics (min_dist)
#The greater min_dist the less confidence we have in the model’s prediction.
#DID NOT WORK; look into it later (AT), the value is in nn_stats
#sc.pl.umap(adams, color='min_dist', vmax=.1, save=imag_name4)

#we input a query profile representing the cell state of interest, querycell state
fm_basic_signature = ['SPP1', 'TREM2', 'GPNMB', 'MMP9', 'CHIT1', 'CHI3L1']
sc.tl.score_genes(adams, fm_basic_signature)
#sc.pl.umap(adams, color='score', save=imag_name5)      (uncomment LATER (AT))

"""
######################################## QUERY MODEL #####################################
"""

# Select the top scoring cells to define our query cell state
sig_query_threshold = adams.obs.score.quantile(.999)
cells_used_in_query = adams.obs.score>=sig_query_threshold
adams.obs['used_in_query'] = cells_used_in_query

#Compute centroid of top scoring cells and compute embedding
from scimilarity.utils import get_centroid
avg_cell = get_centroid(adams.layers['counts'][adams.obs['used_in_query']])
avg_embedding = ca.get_embeddings(avg_cell)

from scimilarity import CellQuery

query_path = '../../models/query_model_v1'
cq = CellQuery(model_path=annotation_path,
               cellsearch_path=query_path)

k = 10000
nn_idxs, nn_dists, results_metadata = cq.search(avg_embedding, k=k)

def calculate_disease_proportions(metadata):
    study_proportions = metadata.disease.value_counts()
    return 100*study_proportions / study_proportions.sum()

def plot_proportions(df, title=None):
    ax = df.plot(kind='barh',
            xlabel='percent of cells',
            title=title,
            grid=False,
            figsize=(4,4))
    ax.tick_params(axis='y', labelsize=8)
    ax.set_xticklabels([f'{int(tick)}%' for tick in ax.get_xticks()]);

query_study = 'DS000011735'
filtered_result_metadata = results_metadata[results_metadata.study!=query_study]
query_disease_frequencies = calculate_disease_proportions(filtered_result_metadata)
plot_proportions(query_disease_frequencies,
                 title='disease proportions for most similar cells')

print("done")
