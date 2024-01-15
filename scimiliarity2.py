""" this is a redundant file but makes a figure with 0 and 1"""

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

# Select visualization options
ALL = 0
clstr = 0
uncnstrnd_prdctn = 0
cnstrnd_prdctn = 0
min_dist = 1
gene_score = 1

annotation_path = '../models/annotation_model_v1'
#data_path2 = '../data/GSE136831_subsample.h5ad'
data_path2 = '../data/KPTracer-Data/expression/adata_processed_combined_SCANTN.h5ad' #(AT)
data_path = '../data/KPTracer-Data/expression/adata_processed_comp_SCANTN.h5ad' #(AT)
imag_name1 = utils_AT.create_image_name("_KPTrcrClstrCOMP_")
imag_name2 = utils_AT.create_image_name("_KPTrcrConstrndAntnCOMP_")
imag_name3 = utils_AT.create_image_name("_KPTrcrConstrneAntnCOMP_")
imag_name4 = utils_AT.create_image_name("_KPTrcrMinDistCOMP_")
imag_name5 = utils_AT.create_image_name("_KPTrcrGeneScoreCOMP_")
imag_name6 = utils_AT.create_image_name("_KPTrcrGeneScoreDscrtCOMP_")

""""  ###### VARIABLES ########  """

# Model read, data read
sc.set_figure_params(dpi_save=800)
warnings.filterwarnings('ignore')
#ca = CellAnnotation(model_path=annotation_path)
# adams = sc.read(data_path)
# adams2 = sc.read(data_path2)

# adams = sc.read(data_path2)
adams_comp = sc.read(data_path)

""" heatmap in yang vs scimilarity """
import seaborn as sns
import matplotlib.pyplot as plt

column1_data = adams_comp.obs['Cluster-Name']
column2_data = adams_comp.obs['celltype_hint']

# Create a frequency matrix
frequency_matrix = pd.crosstab(column1_data, column2_data)

# Create a custom annotation matrix where all values are displayed as strings
annotation_matrix = frequency_matrix.astype(str)

# Create a heatmap using seaborn
# sns.heatmap(frequency_matrix, cmap='coolwarm', annot=True, fmt="d", linewidths=0.5, cbar=False)
sns.heatmap(frequency_matrix, cmap='coolwarm', annot=annotation_matrix, fmt="", linewidths=0.5, cbar=False)

plt.title('Heatmap Between Column1 and Column2')
plt.savefig('heatmap.png', bbox_inches='tight', dpi=150)
plt.show()

print()









fm_basic_signature = (utils_AT.read_csv_column
                      ("../files/yangCOMP_TargetGene.csv", "Common_gene_name"))
sc.tl.score_genes(adams_comp, fm_basic_signature)
if ALL or gene_score:
    sc.pl.umap(adams_comp, color='score', save=imag_name5)

# make a subset of anndata which gene score more than 0.5
threshold = 0.5
values = adams.obs['score'].values
colors = np.where(values > threshold, 'blue', 'gray')
adams.obs['gene_score_dscrt'] = colors

# we are creating a subset just with score more than threlhold
subset_ftr = ['blue']
adamsUPgnscr = adams[adams.obs.iloc[:, -1].isin(subset_ftr)]


""" print out unconstrained predictions in a csv file """
# celltype_counts = adams.obs.predictions_unconstrained.value_counts()
# max_cellcount = celltype_counts
# import csv
# csv_filename = ("../files/yang_uncnstrnd_prdctn.csv")
# cell_name_arr = []
# cell_num_arr = []
# for i in range(max_cellcount.shape[0]):
#         #final_data = (max_cellcount.index[i], max_cellcount[i])
#         cell_name_arr.append(max_cellcount.index[i])
#         cell_num_arr.append(max_cellcount[i])
# with open(csv_filename, 'w', newline='') as csvfile:
#     csv_writer = csv.writer(csvfile)
#     # Write header
#     csv_writer.writerow(["YANG_data_scimilarity_unconstrained_prediction", "Numbers"])
#     for cellName, cellNum in zip(cell_name_arr, cell_num_arr) :
#          csv_writer.writerow([cellName, cellNum])
#     #csv_writer.writerow(zip(cell_name_arr, cell_number_arr))
# print()

""" plot the cell annotation distribution """
import matplotlib.pyplot as plt

# Load the data from the CSV file
df = pd.read_csv("../files/yang_uncnstrnd_prdctn.csv")

fig, ax = plt.subplots()
ax.barh(df['YANG_data_scimilarity_unconstrained_prediction'][:10], df['Ratio'][:10])
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel("Cell Ratio (%)")
ax.set_title("Yang data scimilarity annotation distribution")
image_name = "figures/"+utils_AT.create_image_name("Yang_dist_ratio", format='.jpg')
plt.savefig(image_name, bbox_inches='tight', dpi=150)
plt.show()
print()



import matplotlib.pyplot as plt
import seaborn as sns

x = adams_comp.obs["min_dist"]
x = pd.DataFrame(x)
x["data_type"]="comprehensive"

y = adams.obs["min_dist"]
y = pd.DataFrame(y)
y["data_type"]="NOT_comprehensive"

x = pd.concat([x, y], axis=0)

plt.figure(figsize=(8,6))
sns.violinplot(x = "data_type", y="min_dist", data=x)
plt.grid()
plt.xlabel('Data_type')
plt.ylabel('Min_dist')
plt.title('Scimilarity min distance for Yang datatypes')
plt.savefig(('figures/'+utils_AT.create_image_name('min_dist_distribution_',format='.jpg')), bbox_inches='tight', dpi=150)
plt.show()



sc.pl.umap(adams_comp, color='Cluster-Name', legend_fontsize=4, save=imag_name1, legend_loc='on data')
imag_name1 = utils_AT.create_image_name("_KPTrcrClstrCOMP_")
sc.pl.umap(adams_comp, color='Cluster-Name', legend_fontsize=5, save=imag_name1)
imag_name1 = utils_AT.create_image_name("_KPTrcrClstrCOMP_")
sc.pl.umap(adams_comp, color='celltype_hint', legend_fontsize=5, save=imag_name1)


import matplotlib.pyplot as plt

""" Gene Score Signature """

# we input a query profile representing the cell state of interest, querycell state
#fm_basic_signature = ['SPP1', 'TREM2', 'GPNMB', 'MMP9', 'CHIT1', 'CHI3L1']
fm_basic_signature = (utils_AT.read_csv_column
                      ("../files/common_adata_targetGene_2023_9_20_20_54_19.csv", "Common_gene_name"))
sc.tl.score_genes(adams, fm_basic_signature)
# if ALL or gene_score:
#     sc.pl.umap(adams, color='score', save=imag_name5)

# make a subset of anndata which gene score more than 0.5
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


adams.obs['SC_difrnt_clstr'] = adams_comp.obs['celltype_hint'].values
cell_loc = np.where(adams.obs['SC_difrnt_clstr'].values.codes != adams.obs['celltype_hint'].values.codes, 0, 1)

umap1 = adams.obsm['X_umap'][:, 0]
umap2 = adams.obsm['X_umap'][:, 1]
color_map = {1: 'gray', 0: 'blue'}
colors = [color_map[val] for val in cell_loc]
plt.scatter(umap1, umap2, c=colors)
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.title('Different annotation after comprehensive gene info')
plt.savefig("figures/different_cell_antn_old2.png")
plt.show()
print()