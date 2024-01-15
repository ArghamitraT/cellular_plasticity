# Environment settings
import scanpy as sc
import utils_AT

sc.set_figure_params(dpi=500)
import warnings
warnings.filterwarnings('ignore')
from scimilarity.utils import lognorm_counts
from scimilarity import CellAnnotation, align_dataset


# Instantiate the CellAnnotation object.
annotation_path = '../../models/annotation_model_v1'
ca = CellAnnotation(model_path=annotation_path)


#data_path = '../data/GSE136831_subsample.h5ad'
data_path = '../../data/KPTracer-Data/expression/adata_processed.combined.h5ad'  #(AT)
adams = sc.read(data_path)

#(AT) to match the target gene, all name of the gene needs to be uppercase
for indx_num in range(len(adams.var.index)):
     adams.var.index.values[indx_num] = adams.var.index[indx_num].upper()


adams = align_dataset(adams, ca.gene_order)
adams = lognorm_counts(adams)
adams.obsm['X_scimilarity'] = ca.get_embeddings(adams.X)

#(AT)
image_name = utils_AT.create_image_name("_Dyang3732Lkb1T1_")
sc.pp.neighbors(adams, use_rep='X_scimilarity')

#only show the cells with tumor 3732_Lkb1_T1
subset_ftr = ['3732_Lkb1_T1']
subset_adams = adams[adams.obs.Tumor.isin(subset_ftr)]

sc.tl.umap(subset_adams)

#sc.pl.umap(adams, color='celltype_raw', legend_fontsize=5)
sc.pl.umap(adams, color='Cluster-Name', legend_fontsize=5, save=image_name, legend_loc = "on data") #(AT)
print("done")