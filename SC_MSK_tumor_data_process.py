import scanpy as sc
def split_anndata(adata, chunk_size):
    n_cells = adata.shape[0]
    chunks = []

    for start_idx in range(0, n_cells, chunk_size):
        end_idx = min(start_idx + chunk_size, n_cells)
        chunk = adata[start_idx:end_idx].copy()
        chunks.append(chunk)

    return chunks

def save_chunks(chunks, base_filename):
    for i, chunk in enumerate(chunks):
        chunk_filename = f"{base_filename}_chunk{i}.h5ad"
        chunk.write(chunk_filename)
        print(f"Saved chunk {i} to {chunk_filename}")


import scipy.sparse as sp


# Load your dataset
adata = sc.read('../data/MSK_tumor_data/adata.combined.mnnc.010920.h5ad')
adata_X_csr = sp.csr_matrix(adata.X)

# Now 'adata_X_csr' is the CSR matrix version of 'adata.X'
# To replace the original ndarray in 'adata' with this CSR matrix:
adata.X = adata_X_csr

adata.layers['counts'] = adata.X
# Split the AnnData object into chunks
chunk_size = 10000  # You can adjust this number based on your memory capacity
chunks = split_anndata(adata, chunk_size)
base_filename = '../data/MSK_tumor_data/chunk_data/adata_combined_mnnc_010920'
save_chunks(chunks, base_filename)