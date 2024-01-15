# Given all the overlap.csv file this code makes a side by side violin plot for all cell cluster

import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
import pandas as pd
import plasticity_simlib as psl
import utils_AT as util

""""  ###### VARIABLES ########  """

overlap_dir = "../files/overlaps"
gini = 1
avgOvrlp = 0

fixed_distance = 0
fixed_cellNum = 1

""" KP tumor fixed cell number """      # change
all_files = ["overlapfixedcellnum_NEW_yang_Mesenchymal_1_3457_Apc_2023_11_8_21_24_35.csv",
"overlapfixedcellnum_NEW_yang_Early_gastric_3457_Apc_2023_11_8_21_24_35.csv"]

""" KP tumor fixed distance """         # change
# all_files = ["overlap_NEW_yang_Mesenchymal_1_3457_Apc_2023_11_8_21_17_59.csv",
#  "overlap_NEW_yang_Early_gastric_3457_Apc_2023_11_8_21_17_59.csv"]

""" Other cell cluster """         # change
# all_files = ["overlap_NEW_yang_GASTRIC_secretoryCell_2023_10_15_18_22_37.csv",
#              "overlap_NEW_yang_HP_scrtrCell_2023_10_15_18_33_28.csv",
#              "overlap_NEW_yang_PrEMT_scrtrCell_2023_10_15_18_33_28.csv",
#              "overlap_NEW_yang_Mesen2_cardiacCell_2023_10_15_18_42_15.csv",
#              "overlap_NEW_yang_Mesen2Met_cardiacCell_2023_10_15_18_42_15.csv",
#              "overlap_NEW_yang_Mesen2_epthCell_2023_10_15_18_50_51.csv",
#              "overlap_NEW_yang_Mesen2Met_epthCell_2023_10_15_18_50_51.csv",
#              "overlap_NEW_yang_PrEMT_glndEpthCell_2023_10_16_11_04_35.csv",
#              "overlap_NEW_yang_HP_glndEpthCell_2023_10_16_11_04_35.csv"]

""""  ###### VARIABLES ########  """


i = 0
for file in all_files:
    overlaps = pd.read_csv(os.path.join(overlap_dir, file))
    overlaps = overlaps.set_index(overlaps.columns[0])
    overlaps.columns = overlaps.columns.astype(int)

    # if fixed_cellNum:
    #     """ neighborhood overlap (fixed cell number) """ # change
    #     data = overlaps.melt()
    #     data['value'] = data['value'] / data.iloc[-1, 1]
    #     sns.lineplot(x='variable', y='value', data=data, errorbar='sd')
    #     plt.xlabel('Radius')
    #     plt.ylabel('Neighbourhood Overlap (fixed cell number)')
    #     plt.grid()
    #     name_str = (file.split(".")[0]).split("_")
    #     title =  name_str[3] + "_" + name_str[4]
    #     plt.title(title)
    #     plt.savefig(('figures/' + util.create_image_name('neighborOverlap_' + title, format='.jpg')), bbox_inches='tight',
    #                 dpi=150)
    #     plt.show()
    #
    # if fixed_distance:
    #     """ neighborhood overlap (fixed distance) """  # change
    #     data = overlaps.melt()
    #     data['value'] = data['value'] / data.iloc[-1, 1]
    #     sns.lineplot(x='variable', y='value', data=data, errorbar='sd')
    #     plt.xlabel('Radius')
    #     plt.ylabel('Neighbourhood Overlap (fixed distance)')
    #     plt.grid()
    #     name_str = (file.split(".")[0]).split("_")
    #     title = name_str[3] + "_" + name_str[4]
    #     plt.title(title)
    #     plt.savefig(('figures/' + util.create_image_name('neighborOverlap_' + title, format='.jpg')), bbox_inches='tight',
    #                 dpi=150)
    #     plt.show()

    if avgOvrlp:
        """ (1-avg overlap) """ # change
        overlaps = overlaps / overlaps.shape[0]
        cell_names = overlaps.index
        overlaps = overlaps.values
        gini_I = np.array([1-np.mean(overlaps[i]) for i in range(overlaps.shape[0])])
        temp_gini_index = pd.DataFrame(gini_I, index=cell_names, columns=['Gini Index'])

    if gini:
        """ gini index """ # change
        temp_gini_index = (psl.get_gini_index(overlaps))

    # naming
    name_str = (file.split(".")[0]).split("_")
    temp_gini_index['Cell_type'] = name_str[3]+"_"+name_str[4]+"_"+name_str[5]+"_"+name_str[6]
    if i==0:
        gini_index = temp_gini_index
    else:
        gini_index = pd.concat([gini_index, temp_gini_index], axis=0)
    i+=1
    print()

if fixed_cellNum:
    if gini:
        """ plot Gini index (fixed cell number) """ # change
        plt.figure(figsize=(8,6))
        sns.violinplot(x='Cell_type', y='Gini Index', data=gini_index)
        plt.xticks(rotation=90)
        plt.xlabel('Cell Cluster')
        plt.ylabel('Gini-index')
        plt.grid()
        plt.title("Gini index of different cells (fixed cell number)")
        plt.savefig(('figures/'+util.create_image_name('GI_KPtumor',format='.png')), bbox_inches='tight', dpi=150)
        plt.show()

    if avgOvrlp:
        """ plot 1-avg overlap (fixed cell number) """  # change
        plt.figure(figsize=(8,6))
        sns.violinplot(x='Cell_type', y='Gini Index', data=gini_index)
        plt.xticks(rotation=90)
        plt.xlabel('Cell Cluster')
        plt.ylabel('1-(avg overlap)~plasticity')
        plt.grid()
        plt.title("Gini-ish index of different cells (fixed cell number)")
        plt.savefig(('figures/'+util.create_image_name('GI_ish_KPtumor',format='.png')), bbox_inches='tight', dpi=150)
        plt.show()

if fixed_distance:
    if gini:
        """ plot Gini index (fixed distance) """ # change
        plt.figure(figsize=(8,6))
        sns.violinplot(x='Cell_type', y='Gini Index', data=gini_index)
        plt.xticks(rotation=90)
        plt.xlabel('Cell Cluster')
        plt.ylabel('Gini-index')
        plt.grid()
        plt.title("Gini index of different cells (fixed distance)")
        plt.savefig(('figures/'+util.create_image_name('GI_KPtumor',format='.png')), bbox_inches='tight', dpi=150)
        plt.show()

    if avgOvrlp:
        """ plot 1-avg overlap (fixed distance) """ # change
        plt.figure(figsize=(8,6))
        sns.violinplot(x='Cell_type', y='Gini Index', data=gini_index)
        plt.xticks(rotation=90)
        plt.xlabel('Cell Cluster')
        plt.ylabel('1-(avg overlap)~plasticity')
        plt.grid()
        plt.title("Gini-ish index of different cells (fixed distance)")
        plt.savefig(('figures/'+util.create_image_name('GI_ish_KPtumor',format='.png')), bbox_inches='tight', dpi=150)
        plt.show()
print()
