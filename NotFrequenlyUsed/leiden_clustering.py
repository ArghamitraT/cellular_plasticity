import scanpy as sc
import utils_AT

# Load your data (replace 'your_file.h5ad' with your data file)
adata = sc.read('../../data/KPTracer-Data/expression/adata_processed.combined.h5ad')

sc.pp.filter_cells(adata, min_genes=200)  # Filter out cells with too few genes
sc.pp.filter_genes(adata, min_cells=3)    # Filter out genes expressed in too few cells
sc.pp.normalize_total(adata)              # Normalize total counts per cell
sc.pp.log1p(adata)                        # Log-transform the data

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)  # Compute a neighborhood graph
sc.tl.leiden(adata, resolution=0.5)

img_name = utils_AT.create_image_name("_leiden_")
sc.tl.umap(adata)  # Compute UMAP for visualization
sc.pl.umap(adata, color=['leiden'], save=img_name)