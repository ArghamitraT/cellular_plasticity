# VAE step 2: find the query cell's 50 nearest neighbours in the SCimilarity reference.
#
# Uses the query cell picked in step 1 (results/vae_steps/current_query_cell.txt);
# reads and writes that cell's folder, results/vae_steps/<query cell>/.
#
# WORKFLOW
#   Step 1. Load the query cell's SCimilarity embedding (128 numbers), computed earlier by
#           annotation_entropy_vary_k.py (raw counts, normalized once).
#   Step 2. Search the SCimilarity reference index for its 50 nearest reference cells
#           (cosine distance, same search SCimilarity uses to annotate cells).
#   Step 3. Print every neighbour (reference cell ID, cell type, distance), then the %
#           of each cell type.
#   Step 4. Check the % match step 1 (which read the saved neighbours from the earlier run).
#
# Output (results/vae_steps/<query cell>/):
#   2_neighbours.csv         the 50 neighbours, closest first
#   2_neighbour_percent.csv  % of each cell type among the 50
#
# Run: python 2_query_cell_neighbours.py   (~2 min, mostly loading the 9.3 GB index)

import os

# Numba's default OpenMP threading layer segfaults in this env (pip numba vs
# conda llvm-openmp). scimilarity imports scanpy, so set this before importing it.
os.environ["NUMBA_THREADING_LAYER"] = "workqueue"

import numpy as np
import pandas as pd
from scimilarity import CellAnnotation

pd.set_option("display.width", 200)


# --------------------------------------------------
# Settings
# --------------------------------------------------

root = "/Users/shivanim/Documents/peer_lab/scimilarity_test"
model_path = os.path.join(root, "models/model_v1.1")
nn_dir = os.path.join(root, "results/annotation_entropy_vary_k")  # from annotation_entropy_vary_k.py
steps_dir = os.path.join(root, "results/vae_steps")

k = 50


# --------------------------------------------------
# Step 1. The query cell's embedding
# --------------------------------------------------

with open(os.path.join(steps_dir, "current_query_cell.txt")) as f:
    query_cell = f.read().strip()
out_dir = os.path.join(steps_dir, query_cell)
step1 = pd.read_csv(os.path.join(out_dir, "1_query_cell.csv"))

cells = pd.read_csv(os.path.join(nn_dir, "cells_all_true_counts.csv"), index_col=0)
z_query = np.load(os.path.join(nn_dir, "embeddings_all_true_counts.npy"))[[cells.index.get_loc(query_cell)]]
print(f"Query cell: {query_cell} ({cells.loc[query_cell, 'cell_type']}), "
      f"embedding of {z_query.shape[1]} numbers")


# --------------------------------------------------
# Step 2. Its 50 nearest reference cells
# --------------------------------------------------

ref_labels = pd.read_csv(os.path.join(model_path, "annotation/reference_labels.tsv"),
                         sep="\t", header=None, usecols=[0])[0]

print("Loading SCimilarity reference index (labelled_kNN.bin, 9.3 GB, takes a minute)...")
ca = CellAnnotation(model_path)
_, nn_idxs, nn_dists, _ = ca.get_predictions_knn(z_query, k=k)
del ca  # free the memory

neighbours = pd.DataFrame({"rank": np.arange(1, k + 1),
                           "reference_cell_id": nn_idxs[0],
                           "cell_type": ref_labels.values[nn_idxs[0]],
                           "cosine_distance": np.round(nn_dists[0], 4)})
neighbours.to_csv(os.path.join(out_dir, "2_neighbours.csv"), index=False)


# --------------------------------------------------
# Step 3. Print neighbours and %
# --------------------------------------------------

print(f"\nIts {k} nearest SCimilarity reference cells (closest first):")
print(neighbours.to_string(index=False))

percent = neighbours["cell_type"].value_counts().rename_axis("cell_type").reset_index(name="n_neighbours")
percent["pct_of_neighbours"] = 100 * percent["n_neighbours"] / k
percent.to_csv(os.path.join(out_dir, "2_neighbour_percent.csv"), index=False)

print(f"\n% of each cell type among the {k} neighbours:")
for _, r in percent.iterrows():
    print(f"  {r['pct_of_neighbours']:5.1f}%  ({r['n_neighbours']:2d} of {k})  {r['cell_type']}")


# --------------------------------------------------
# Step 4. Same as step 1?
# --------------------------------------------------

before = step1.set_index("type")["n_neighbours"]
now = percent.set_index("cell_type")["n_neighbours"]
same = before.sort_index().equals(now.sort_index())
print(f"\nSame neighbour counts as step 1: {same}")
if not same:
    print(pd.DataFrame({"step1": before, "now": now}).fillna(0).astype(int))
print("Saved in", out_dir)
