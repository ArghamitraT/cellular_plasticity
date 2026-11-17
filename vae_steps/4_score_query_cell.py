# VAE step 4: score the query cell on each neighbour type's VAE.
#
# Uses: the query cell (step 1, results/vae_steps/current_query_cell.txt), its neighbour %
# (step 2), the shared trained VAEs and held-out real cells (step 3,
# results/vae_steps/3_vaes_<run_name>/). Nothing is trained here.
#
# HOW THE SCORE IS CALCULATED: score() of the VAE version chosen in vae_model.py
#   "learned_variance": log-likelihood, estimated by importance sampling
#                       (vae_model_learned_variance.py)
#   "mse":              -(MSE reconstruction + KL), with h = mu (vae_model_mse.py)
#   Higher = the VAE explains the cell well = looks like a normal cell of the type.
#   Scores can only be compared with other scores from the same VAE.
#
# WORKFLOW
#   Step 1. Load the query cell's embedding, its neighbour %, the VAEs and held-out cells.
#   Step 2. Score the held-out real cells of each type under that type's VAE
#           = what a normal cell of that type scores (grey histogram).
#   Step 3. Score the query cell under each type's VAE.
#   Step 4. Table, per type (each computed from that type's VAE alone):
#           - fit percentile: % of held-out real cells the query scores higher than
#             (50% = typical, ~0% = worse than nearly all real cells)
#           - rank of each type by neighbours and by VAE score
#   Step 5. Histograms: one panel per type, grey = real cells, red line = query cell.
#
# Output (results/vae_steps/<query cell>/<run_name>/):
#   4_results.csv, 4_histograms.png,
#   run_info.txt (which VAE version and settings made the results; a line added per run)
#
# Run: python 4_score_query_cell.py   (under a minute)

import os

import numpy as np
import pandas as pd
import torch
from matplotlib import pyplot as plt

from vae_model import VAE, run_name, score_label, write_run_info

pd.set_option("display.width", 200)


# --------------------------------------------------
# Settings
# --------------------------------------------------

root = "/Users/shivanim/Documents/peer_lab/scimilarity_test"
nn_dir = os.path.join(root, "results/annotation_entropy_vary_k")  # query embeddings
steps_dir = os.path.join(root, "results/vae_steps")
vae_dir = os.path.join(steps_dir, f"3_vaes_{run_name}")  # one folder per VAE version


def file_name(t):
    return t.replace(" ", "_").replace(",", "")


def pct_text(p):
    return "<0.1%" if p < 0.1 else f"{p:.1f}%"


# --------------------------------------------------
# Step 1. Load everything
# --------------------------------------------------

with open(os.path.join(steps_dir, "current_query_cell.txt")) as f:
    query_cell = f.read().strip()
cell_dir = os.path.join(steps_dir, query_cell)   # steps 1-2 results (read)
out_dir = os.path.join(cell_dir, run_name)       # steps 3-4 results (written)
os.makedirs(out_dir, exist_ok=True)
cells = pd.read_csv(os.path.join(nn_dir, "cells_all_true_counts.csv"), index_col=0)
z_query = np.load(os.path.join(nn_dir, "embeddings_all_true_counts.npy"))[[cells.index.get_loc(query_cell)]]

pct = pd.read_csv(os.path.join(cell_dir, "2_neighbour_percent.csv"))
types = [t for t in pct["cell_type"] if os.path.exists(os.path.join(vae_dir, file_name(t) + ".pt"))]
print(f"Query cell: {query_cell}")
print(f"VAE version: {run_name}")
print(f"Scoring it on {len(types)} VAEs: {', '.join(types)}")

vaes, x_held = {}, {}
for t in types:
    vaes[t] = VAE()
    vaes[t].load_state_dict(torch.load(os.path.join(vae_dir, file_name(t) + ".pt")))
    vaes[t].eval()
    x_held[t] = np.load(os.path.join(vae_dir, file_name(t) + "_heldout.npy"))


# --------------------------------------------------
# Steps 2-3. Score real cells and the query cell
# --------------------------------------------------

held_score, query_score = {}, {}
for t in types:
    held_score[t] = vaes[t].score(x_held[t])
    query_score[t] = float(vaes[t].score(z_query)[0])


# --------------------------------------------------
# Step 4. Table
# --------------------------------------------------

res = pct.set_index("cell_type").loc[types, ["n_neighbours", "pct_of_neighbours"]].reset_index()
res = res.rename(columns={"cell_type": "type"})
res["query_score"] = [query_score[t] for t in types]
res["typical_real_cell_score"] = [float(np.median(held_score[t])) for t in types]
res["fit_percentile"] = [100 * (held_score[t] < query_score[t]).mean() for t in types]
res["rank_by_neighbours"] = res["n_neighbours"].rank(ascending=False, method="min").astype(int)
res["rank_by_vae"] = res["query_score"].rank(ascending=False, method="min").astype(int)
res = res.round(4)  # scores are small (around -0.1), so keep 4 decimals
res.insert(0, "vae_version", run_name)
res.to_csv(os.path.join(out_dir, "4_results.csv"), index=False)

show = pd.DataFrame({
    "type": res["type"],
    "neighbours": [f"{p:.0f}%" for p in res["pct_of_neighbours"]],
    "query score": res["query_score"],
    "typical real cell": res["typical_real_cell_score"],
    "query beats": [pct_text(p) + " of real cells" for p in res["fit_percentile"]],
    "rank (neighbours)": res["rank_by_neighbours"],
    "rank (VAE)": res["rank_by_vae"],
})
print(f"\nResults (score = {score_label}, higher = fits the type better):")
print(show.to_string(index=False))


# --------------------------------------------------
# Step 5. Histograms
# --------------------------------------------------

# Panels stacked top to bottom, all on the same x-axis so positions compare directly.
# x-axis: from the lowest to the highest score, plus a 5% margin on each side.
all_scores = np.concatenate([held_score[t] for t in types] + [res["query_score"].values])
margin = 0.05 * (all_scores.max() - all_scores.min())
x_lo, x_hi = all_scores.min() - margin, all_scores.max() + margin
bins = np.linspace(x_lo, x_hi, 300)  # same bins in every panel

fig, axes = plt.subplots(len(types), 1, figsize=(12, 2.4 * len(types)), sharex=True)
for ax, (_, r) in zip(np.atleast_1d(axes), res.iterrows()):
    t = r["type"]
    ax.hist(held_score[t], bins=bins, color="grey",
            label=f"real cells of this type (n = {len(held_score[t])})")
    ax.axvline(r["typical_real_cell_score"], color="black", linestyle="--", linewidth=1,
               label="typical real cell")
    ax.axvline(r["query_score"], color="red", linewidth=2.5, label="SCLC query cell")
    ax.set_title(f"{t}  ({r['pct_of_neighbours']:.0f}% of neighbours)", fontsize=11, loc="left")
    ax.text(0.01, 0.92, f"grey: {len(held_score[t])} real cells\n"
                        f"query: {r['query_score']:.4f}\n"
                        f"typical real cell: {r['typical_real_cell_score']:.4f}\n"
                        f"query beats {pct_text(r['fit_percentile'])} of real cells",
            transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(facecolor="white", edgecolor="lightgrey"))
    ax.set_ylabel("cells")
    ax.grid(axis="x", alpha=0.3)
np.atleast_1d(axes)[0].legend(loc="upper center", fontsize=8)
np.atleast_1d(axes)[-1].set_xlim(x_lo, x_hi)
np.atleast_1d(axes)[-1].set_xlabel(f"{score_label}  (higher = fits the type better)")
fig.suptitle(f"SCLC-A query cell {query_cell}: fit to each neighbour type's VAE ({run_name})",
             fontsize=13)
fig.text(0.5, -0.01,
         "Grey = real reference cells of that type that the VAE never trained on (n in each panel). "
         "Red line = the SCLC query cell.\nFurther left = fits worse. "
         "Panels ordered by % of the query cell's 50 SCimilarity neighbours; all share the same x-axis.",
         ha="center", va="top", fontsize=10)
fig.tight_layout()
fig_path = os.path.join(out_dir, "4_histograms.png")
fig.savefig(fig_path, bbox_inches="tight", dpi=150)
plt.close(fig)
write_run_info(out_dir, query_cell, "4_score_query_cell.py")
print("\nSaved:", os.path.join(out_dir, "4_results.csv"))
print("Saved:", fig_path)
print("Which VAE version made these results:", os.path.join(out_dir, "run_info.txt"))
