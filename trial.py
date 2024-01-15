# Environment settings
import scanpy as sc
import warnings
import os
import numpy as np

#data_path = '../data/MSK_tumor_data/glasner_etal_globalAnndata_20230112.vHTA.h5ad'
data_path = '../data/MSK_tumor_data/adata.combined.mnnc.010920.h5ad'
adams = sc.read(data_path)


# Replace 'file_path.npy' with the path to your .npy file
data = np.load('../models/query_model_v1/test_embedding.npy')

# Now `data` holds the NumPy array that was stored in the .npy file
print(data)


# Define directory names
dir1 = "path/to/first_directory"
dir2 = "path/to/second_directory"

# Concatenate directory names
concatenated_dir = os.path.join(dir1, dir2)

# Print the concatenated directory
print(concatenated_dir)


data_path = '../data/KPTracer-Data/expression/adata_processed.combined.h5ad' #(AT)

sc.set_figure_params(dpi=100)
warnings.filterwarnings('ignore')
adams = sc.read(data_path)

print("done")










# # Environment settings
# import numpy as np
# import scanpy as sc
# #sc.set_figure_params(dpi=100)
#
# # import warnings
# # warnings.filterwarnings('ignore')
# #
# # from scimilarity.utils import lognorm_counts
# # from scimilarity import CellAnnotation, align_dataset
# #
# # # Instantiate the CellAnnotation object.
# # annotation_path = '../models/annotation_model_v1'
# # ca = CellAnnotation(model_path=annotation_path)
#
# # Load the tutorial data.
# data_path = ('../data/KPTracer-Data/expression/adata_processed.combined.h5ad')
# adams_dyang = sc.read(data_path)
#
# data_path = '../data/GSE136831_subsample.h5ad'
# adams_atlas = sc.read(data_path)
#
# vars_dyang = []
# vars_atlas = []
# vars_dyang = adams_dyang.var_names
# vars_atlas = adams_atlas.var_names
#
# for dyang in vars_dyang:
#     for atlas in vars_atlas:
#         if (dyang == atlas):
#             print("found ", dyang)
#
#
# print("done")