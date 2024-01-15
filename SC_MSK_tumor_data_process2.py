""" We are processing the chunked Tabula Sapiens data """

import os
import scanpy as sc
from scimilarity.utils import lognorm_counts
from scimilarity import CellAnnotation, align_dataset

annotation_path = '../models/annotation_model_v1'
ca = CellAnnotation(model_path=annotation_path)

# Define a function to process chunks
def process_chunk(adams):
    adams = align_dataset(adams, ca.gene_order)
    adams = lognorm_counts(adams)
    adams.obsm['X_scimilarity'] = ca.get_embeddings(adams.X)

    # Clustering
    # sc.pp.neighbors(adams, use_rep='X_scimilarity')
    # sc.tl.umap(adams)
    return adams

# Split the AnnData object into chunks

def load_process_save_chunks(chunks_dir, base_filename, chunk_size, save_processed=True):
    processed_chunks = []
    n_chunks = len([name for name in os.listdir(chunks_dir) if
                    os.path.isfile(os.path.join(chunks_dir, name)) and name.startswith(base_filename) and name.endswith(
                        '.h5ad')])

    for i in range(n_chunks):
        chunk_filename = os.path.join(chunks_dir, f"{base_filename}_chunk{i}.h5ad")
        chunk = sc.read(chunk_filename)

        processed_chunk = process_chunk(chunk)
        processed_chunks.append(processed_chunk)
        print("processing done chunk ", i)

    concatenated_adata = sc.concat(processed_chunks, join='outer')
    (concatenated_adata.write
     ('../data/MSK_tumor_data/adata_combined_mnnc_010920_Scimilarity.h5ad'))

        # if save_processed:
        #     processed_chunk.write(os.path.join(chunks_dir, f"{base_filename}_processed_chunk{i}.h5ad"))
        # else:
        #     processed_chunks.append(processed_chunk)

    #return processed_chunks


# Example usage
chunks_dir = '../data/MSK_tumor_data/chunk_data/'
base_filename = 'adata_combined_mnnc_010920'
chunk_size = 10000
# processed_chunks = load_process_save_chunks(chunks_dir, base_filename, chunk_size, save_processed=True)
load_process_save_chunks(chunks_dir, base_filename, chunk_size, save_processed=True)


# Assuming 'chunks' is a list of AnnData objects split from 'adata'
# processed_chunks = [process_chunk(chunk) for chunk in chunks]
#
# # Combine processed chunks back into one dataset
# adata_processed = pd.concat(processed_chunks, axis=0)
# adata.write('../data/TS_filteredCellGene.h5ad')

# sc.settings.set_figure_params(dpi=80, facecolor='white')
# adata = sc.read("../data/aad23e04-0dfd-4826-b150-4362f3936849.h5ad")
#
# adata.var_names_make_unique()
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
#
# adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
# sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
# sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
#              jitter=0.4, multi_panel=True)
#
# adata = adata[adata.obs.n_genes_by_counts < 2500, :]
# adata = adata[adata.obs.pct_counts_mt < 5, :]
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)
# adata.write('../data/TS_filteredCellGene.h5ad')
#
# sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
# sc.pl.highly_variable_genes(adata)
#
# adata = adata[:, adata.var.highly_variable]
# sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
# sc.pp.scale(adata, max_value=10)
# sc.tl.pca(adata, svd_solver='arpack')
# sc.pl.pca(adata, color='CST3')
