"""
THIS FILE CONTAINS SOME SCRATCH CODES THAT COMES HANDY DEALING/MANUPULATING ANNDATA
"""

""" Find out unique values in a column of anndata """
# column_values = adams.obs['predictions_unconstrained']
# unique_values = column_values.unique()

""" print out top 19 unconstrained predictions in a csv file """
# celltype_counts = adams.obs.predictions_unconstrained.value_counts()
# max_cellcount = celltype_counts[:19]
# import csv
# csv_filename = ("../files/top_19_uncnstrnd_prdctn.csv")
# with open(csv_filename, 'w', newline='') as csvfile:
#     csv_writer = csv.writer(csvfile)
#     # Write header
#     csv_writer.writerow(["Top_19_scimilarity_unconstrained_prediction"])
#     for items in max_cellcount :
#         csv_writer.writerow([items])

""" annotation of a particular sub set """
# # cell annotation for high plasticity cells
# subset_ftr = ['High plasticity']
# subset_adams = adams[adams.obs.iloc[:, -2].isin(subset_ftr)]
#
# # Visualization
# sc.pl.umap(subset_adams[subset_adams.obs.predictions_unconstrained.isin(well_represented_celltypes)],
#            color='predictions_unconstrained',
#            legend_fontsize=5, save=imag_name2, legend_loc = "on data")

""" plot the annotation of every Yang's cluster individaully """
# yang_clstrs = adams.obs['Cluster-Name'].unique()
# for clstr in yang_clstrs:
#     subset_adams = adams[adams.obs['Cluster-Name'].isin([clstr])]
#     img_name_clstr = utils_AT.create_image_name("_"+clstr+"_")
#     sc.pl.umap(subset_adams, color='celltype_hint', legend_fontsize=5, save=img_name_clstr, legend_loc="on data", title="Yang_"+clstr)

""" for each Yang annotation, it gives Scimilarity annotations and their %"""
# clstr_antn_prcntg = []
# yang_clstrs = adams.obs['Cluster-Name'].unique()
# for clstr in yang_clstrs:
#     subset_adams = adams[adams.obs['Cluster-Name'].isin([clstr])]
#     column_values = subset_adams.obs['celltype_hint']
#     unique_values = column_values.unique()
#     for celltypes in unique_values:
#         numOfCells = len(subset_adams[subset_adams.obs['celltype_hint'] == celltypes])
#         ratio = round(100 * (numOfCells / subset_adams.shape[0]), 2)
#         final_data = (clstr, celltypes, numOfCells, ratio)
#         clstr_antn_prcntg.append(final_data)
#         print()
# # Specify the CSV file name
# csv_filename = ("../files/Yang_scimilarity_annotation_ratio.csv")
# with open(csv_filename, 'w', newline='') as csvfile:
#     csv_writer = csv.writer(csvfile)
#     csv_writer.writerow(["Yang_annotation", "Scimilarity_annotation", "Num_of_cells_Sc", "Ratio_Sc"])
#     for items in clstr_antn_prcntg:
#         csv_writer.writerow(items)

""" plot the annotations which have overlap gene score more than 0.5 """
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
# # Visualization
# if ALL or clstr:
#     sc.pl.umap(adams, color='Cluster-Name', legend_fontsize=8, save=imag_name1, legend_loc='on data', title="Yang annotation (overlapping gene score>0.5)")

""" Find out rows without a certain column value"""
# selected_rows = np.where(adams.obs['genotype'] != 'sgNT')[0]
# adams_T = adams[selected_rows, :].copy()

""" Target gene printing: print out the overlapping genes with Scimilarity """
# target_gene_arr = []
    # for tar_gene in target_gene_order:
    #     if tar_gene in data.var.index.values:
    #         target_gene_arr.append(tar_gene)
    #
    # import csv
    # csv_filename = ("../files/yangCOMP_TargetGene.csv")
    # with open(csv_filename, 'w', newline='') as csvfile:
    #     csv_writer = csv.writer(csvfile)
    #     # Write header
    #     csv_writer.writerow(["YANG_COMP_data_TargetGeneList"])
    #     for tar_gene in target_gene_arr:
    #         csv_writer.writerow([tar_gene])

""" makes 0/1 plot """
# adams.obs['SC_difrnt_clstr'] = adams_comp.obs['celltype_hint'].values
# cell_loc = np.where(adams.obs['SC_difrnt_clstr'].values.codes != adams.obs['celltype_hint'].values.codes, 0, 1)
#
# umap1 = adams.obsm['X_umap'][:, 0]
# umap2 = adams.obsm['X_umap'][:, 1]
# color_map = {1: 'gray', 0: 'blue'}
# colors = [color_map[val] for val in cell_loc]
# plt.scatter(umap1, umap2, c=colors)
# plt.xlabel('UMAP1')
# plt.ylabel('UMAP2')
# plt.title('Different annotation after comprehensive gene info')
# plt.savefig("figures/different_cell_antn_old.png")
# plt.show()
# print()

""" bar plot the cell annotation distribution """
# import matplotlib.pyplot as plt
#
# # Load the data from the CSV file
# df = pd.read_csv("../files/yangCOMP_uncnstrnd_prdctn.csv")
# fig, ax = plt.subplots()
# ax.barh(df['YANG_COMP_data_scimilarity_unconstrained_prediction'][:10], df['Ratio'][:10])
# ax.invert_yaxis()  # labels read top-to-bottom
# ax.set_xlabel("Cell Ratio (%)")
# ax.set_title("Yang comprehensive data scimilarity annotation distribution")
# image_name = "figures/"+utils_AT.create_image_name("Yang_COMP_dist_ratio", format='.jpg')
# plt.savefig(image_name, bbox_inches='tight', dpi=150)
# plt.show()
# print()