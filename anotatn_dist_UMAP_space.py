# Environment settings
import scanpy as sc
import utils_AT
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import utils_AT as utils
import pandas as pd
import seaborn as sns

# from bokeh.plotting import figure, show, output_file
# from bokeh.models import HoverTool, ColumnDataSource
# from bokeh.transform import linear_cmap
# from bokeh.palettes import Viridis256


""" variable """
data_path = '../data/MSK_tumor_data/adata_combined_mnnc_010920_SCANTN2.h5ad'
adata = sc.read(data_path)
cell_type = 'pulmonary ionocyte'
centers = [(12.196, 9.617), (17.343, 0.874)]
entropy_d = {}

adata = adata[adata.obs['predictions_unconstrained'] == cell_type].copy() # filter cell type
# sc.pl.umap(adata, color='predictions_unconstrained', legend_fontsize=5, save=utils_AT.create_image_name("_ MSK_tumorPredUncnstrn_"))

""" interactive UMAP """
# umap_coords = adata.obsm['X_umap']
# source = ColumnDataSource(data=dict(
#     x=umap_coords[:, 0],
#     y=umap_coords[:, 1],
#     desc=[f'Cell {i}' for i in range(umap_coords.shape[0])]))
# output_file("umap_pulmonary_ionocyte.html")
# tools = "pan,wheel_zoom,box_select,lasso_select,reset"
# p = figure(tools=tools, title="UMAP Plot", toolbar_location="above")
# p.circle('x', 'y', source=source, size=5)
# hover = HoverTool()
# hover.tooltips = [("Index", "@desc"), ("(x,y)", "($x, $y)")]
# p.add_tools(hover)
# show(p)

def calculate_row_entropy(row):
    row = row[row > 0] # Filter out zero values to avoid log(0)
    return -np.sum(row * np.log(row))

# Loop over each cell type
for idx, (center_x, center_y) in enumerate(centers):

    """ choosing a part of UMAP """
    sc.pl.umap(adata, color='predictions_unconstrained', vmax=.1, show=False)
    umap_coords = adata.obsm['X_umap']
    distances = np.sqrt((umap_coords[:, 0] - center_x) ** 2 + (umap_coords[:, 1] - center_y) ** 2)
    num_cells = 1000
    selected_cells_indices = np.argsort(distances)[:num_cells]
    radius = np.max(distances[selected_cells_indices])
    circle = patches.Circle((center_x, center_y), radius, edgecolor='red', facecolor='none', linewidth=2)
    ax = plt.gca()
    ax.add_patch(circle)
    plt.title('UMAP Plot with Highlighted Region')
    # plt.savefig('figures/' + utils_AT.create_image_name('selected_region'))
    plt.show()

    selected_adata = adata[selected_cells_indices].copy()

    """ entropy and frequency """
    freq_mat = utils.make_frequency_matrix(selected_adata.obs['sc_hits'])
    freq_mat.columns = freq_mat.columns.str.replace(' ', '_')
    row_sums = freq_mat.sum(axis=1)
    normalized_df = freq_mat.div(row_sums, axis=0)
    entropies = normalized_df.apply(calculate_row_entropy, axis=1)
    entropy_list = entropies.tolist()

    type_name = 'x_'+str(center_x)+'_y_'+str(center_y)
    entropy_d[type_name] = entropy_list

""" violin plot """
results = pd.Series(entropy_d)
plt.figure(figsize=(16, 10))
colors = ['skyblue', 'lightgreen', 'violet']  # Define colors for each cluster
sns.violinplot(data=results)
plt.xticks(ticks=range(len(results.index)), labels=results.index)
plt.xlabel("Cell clusters")
plt.ylabel("Annotation entropy")
plt.title("MSK Lung cancer data " + cell_type)
plt.grid()
plt.savefig("figures/"+utils.create_image_name("entropy"))
plt.show()
print()







