# VAE step 3: train one VAE per neighbour cell type of the query cell.
#
# The query cell (results/vae_steps/current_query_cell.txt) and its neighbour types come
# from steps 1-2. Each VAE is trained on SCimilarity REFERENCE cells of its type, not on
# the 50 neighbours: 50 cells are far too few, and the VAE should learn what that type
# normally looks like. Which VAE version is used is set in vae_model.py
# ("learned_variance" or "mse").
#
# VAEs are SHARED between query cells (a VAE depends only on its cell type), in
# results/vae_steps/3_vaes_<run_name>/. A type that already has a VAE is reused, not
# retrained, so different query cells are scored on the same VAEs.
#
# WORKFLOW
#   Step 1. Neighbour types from step 2 (types with >= min_neighbours of the 50).
#   Step 2. Reference embeddings for each new type (128 numbers per cell), from the
#           SCimilarity index. Types with more than max_cells_per_type cells: random sample.
#   Step 3. Split each new type's cells: 80% training, 20% held out. The VAE never sees
#           the held-out cells; in step 4 they show what a real cell of that type scores.
#   Step 4. Train one VAE per new type (model in vae_model.py).
#   Step 5. Training check for all of this query cell's types (new and reused):
#           - mean score (score() of the chosen VAE version) of training cells vs held-out cells. Similar = the
#             VAE learned the cell type; training much higher = it memorised its training cells.
#           - mean reconstruction loss and KL of the held-out cells, separately.
#           - validation curves, 3 rows: reconstruction, KL, and total ELBO = -(recon + KL).
#             Should change, then flatten (= learning has plateaued). Figure starts at
#             step 200: the first drop squashes the y-axis.
#
# Output:
#   results/vae_steps/3_vaes_<run_name>/ (shared)   per type: <type>.pt (trained VAE),
#       <type>_heldout.npy (held-out cells, for step 4), <type>_check.csv, <type>_curve.csv
#   results/vae_steps/reference_embeddings/<type>.npy  (shared by all VAE versions)
#   results/vae_steps/<query cell>/<run_name>/   3_training_check.csv, 3_training_curves.png,
#       run_info.txt (which VAE version and settings made the results; a line added per run)
#
# Run: python 3_train_vaes.py   (~5-10 min per new type; reused types are instant)

import os

# Numba's default OpenMP threading layer segfaults in this env (pip numba vs
# conda llvm-openmp). scimilarity imports scanpy, so set this before importing it.
os.environ["NUMBA_THREADING_LAYER"] = "workqueue"

import numpy as np
import pandas as pd
import torch
from matplotlib import pyplot as plt

from vae_model import recon_label, run_name, train_vae, write_run_info

pd.set_option("display.width", 200)


# --------------------------------------------------
# Settings
# --------------------------------------------------

root = "/Users/shivanim/Documents/peer_lab/scimilarity_test"
model_path = os.path.join(root, "models/model_v1.1")
steps_dir = os.path.join(root, "results/vae_steps")
vae_dir = os.path.join(steps_dir, f"3_vaes_{run_name}")  # one folder per VAE version

min_neighbours = 2          # train VAEs for types with >= this many of the 50 neighbours
max_cells_per_type = 20000  # random sample of big types, keeps training fast
heldout_frac = 0.2

ref_dir = os.path.join(steps_dir, "reference_embeddings")  # same for every VAE version
os.makedirs(ref_dir, exist_ok=True)
os.makedirs(vae_dir, exist_ok=True)
rng = np.random.default_rng(0)
torch.manual_seed(0)

with open(os.path.join(steps_dir, "current_query_cell.txt")) as f:
    query_cell = f.read().strip()
cell_dir = os.path.join(steps_dir, query_cell)   # steps 1-2 results (read)
out_dir = os.path.join(cell_dir, run_name)       # steps 3-4 results (written)
os.makedirs(out_dir, exist_ok=True)
print(f"VAE version: {run_name}")


def path(t, suffix):
    return os.path.join(vae_dir, t.replace(" ", "_").replace(",", "") + suffix)


# --------------------------------------------------
# Step 1. Neighbour types from step 2
# --------------------------------------------------

pct = pd.read_csv(os.path.join(cell_dir, "2_neighbour_percent.csv"))
types = list(pct.loc[pct["n_neighbours"] >= min_neighbours, "cell_type"])
new_types = [t for t in types if not os.path.exists(path(t, ".pt"))]

ref_labels = pd.read_csv(os.path.join(model_path, "annotation/reference_labels.tsv"),
                         sep="\t", header=None, usecols=[0])[0]
ref_counts = ref_labels.value_counts()
print(f"Query cell: {query_cell}")
print("Neighbour types:")
for _, r in pct.iterrows():
    t = r["cell_type"]
    if t not in types:
        status = "skipped: too few neighbours"
    elif t in new_types:
        status = "will train a new VAE"
    else:
        status = "already has a VAE, reusing it"
    print(f"  {r['pct_of_neighbours']:5.1f}%  {t}  ({ref_counts[t]} reference cells)  -> {status}")


# --------------------------------------------------
# Step 2. Reference embeddings for each new type
# --------------------------------------------------

def ref_file(t):
    return os.path.join(ref_dir, t.replace(" ", "_").replace(",", "") + ".npy")


missing = [t for t in new_types if not os.path.exists(ref_file(t))]
if missing:
    from scimilarity import CellAnnotation

    print("\nLoading SCimilarity reference index (labelled_kNN.bin, 9.3 GB, takes a minute)...")
    ca = CellAnnotation(model_path)
    for t in missing:
        ids = np.where(ref_labels.values == t)[0]
        if len(ids) > max_cells_per_type:
            ids = rng.choice(ids, max_cells_per_type, replace=False)
        np.save(ref_file(t), np.array(ca.knn.get_items(ids), dtype=np.float32))
    del ca  # free the memory


# --------------------------------------------------
# Steps 3-4. Split and train each new type
# --------------------------------------------------

if new_types:
    print("\nTraining (a few minutes per type)...")
for t in new_types:
    x = np.load(ref_file(t))
    x = x[rng.permutation(len(x))]
    n_held = int(len(x) * heldout_frac)
    x_held, x_train = x[:n_held], x[n_held:]

    vae, curve = train_vae(x_train, rng)
    torch.save(vae.state_dict(), path(t, ".pt"))
    np.save(path(t, "_heldout.npy"), x_held)
    pd.DataFrame(curve, columns=["step", "val_recon", "val_kl"]).to_csv(path(t, "_curve.csv"), index=False)

    score_train = vae.score(x_train[:n_held]).mean()  # same number of cells as held-out
    score_held = vae.score(x_held).mean()
    with torch.no_grad():
        recon_held, kl_held = vae.loss_per_cell(torch.tensor(x_held), sample=False)
    pd.DataFrame([{"type": t, "reference_cells": ref_counts[t], "cells_used": len(x),
                   "train_cells": len(x_train), "heldout_cells": n_held,
                   "mean_score_train": float(score_train),
                   "mean_score_heldout": float(score_held),
                   "gap": float(score_train - score_held),
                   "heldout_recon": float(recon_held.mean()),
                   "heldout_kl": float(kl_held.mean())}]).to_csv(path(t, "_check.csv"), index=False)
    print(f"  {t}: done")


# --------------------------------------------------
# Step 5. Training check for this query cell's types
# --------------------------------------------------

check = pd.concat([pd.read_csv(path(t, "_check.csv")) for t in types], ignore_index=True)
check.insert(0, "vae_version", run_name)
check.to_csv(os.path.join(out_dir, "3_training_check.csv"), index=False)
print("\nTraining check (score: higher = fits better; gap = train - held-out, small = learned "
      "the type, didn't memorise):")
print(check.round(4).to_string(index=False))

curves = {t: pd.read_csv(path(t, "_curve.csv")).set_index("step") for t in types}
for c in curves.values():
    c["val_elbo"] = -(c["val_recon"] + c["val_kl"])  # ELBO = -(reconstruction + KL)
for part in ["val_recon", "val_kl", "val_elbo"]:
    table = pd.DataFrame({t: c[part] for t, c in curves.items()})
    change = table[table.index % 1000 == 0].diff().iloc[1:]  # flattened = near 0 by the end
    print(f"\nChange in validation {part[4:]} per 1000 steps (should shrink towards 0):")
    print(change.round(4).to_string())

# rows: reconstruction loss, KL, total ELBO; one column per type
rows = [("val_recon", f"{recon_label}\nlower = better"),
        ("val_kl", "KL divergence\nlower = closer to prior"),
        ("val_elbo", "ELBO = -(recon + KL)\nhigher = better")]
fig, axes = plt.subplots(len(rows), len(types), figsize=(4 * len(types), 9.5), squeeze=False)
for j, t in enumerate(types):
    c = curves[t][curves[t].index >= 200]
    for i, (part, label) in enumerate(rows):
        ax = axes[i, j]
        ax.plot(c.index, c[part], marker="o", markersize=2)
        ax.set_title(t if i == 0 else "")
        ax.grid(alpha=0.3)
        if j == 0:
            ax.set_ylabel(label)
        if i == len(rows) - 1:
            ax.set_xlabel("training step")
fig.suptitle(f"VAE validation curves ({run_name}): learning has plateaued when a curve goes flat")
fig.tight_layout()
fig.savefig(os.path.join(out_dir, "3_training_curves.png"), bbox_inches="tight", dpi=150)
plt.close(fig)
write_run_info(out_dir, query_cell, "3_train_vaes.py")
print("\nSaved in", out_dir, "(VAEs in", vae_dir + ")")
print("Which VAE version made these results: run_info.txt in that folder")
