# VAE step 1 of 4: pick an SCLC-A query cell and show its neighbours.
#
# Steps (check each step's output before running the next):
#   1. pick the query cell                          <- this script (~30 s)
#   2. find its 50 SCimilarity neighbours
#   3. train one VAE per neighbour type (shared between query cells)
#   4. score the query cell on each VAE (histograms + table)
#
# Each query cell gets its own folder, results/vae_steps/<query cell>/. This script writes
# the chosen cell to results/vae_steps/current_query_cell.txt; steps 2-4 read it from there.
#
# Two ways to pick (setting `pick`):
#   "representative": entropy near the SCLC-A median (median +- entropy_window), top n_types
#                     neighbour types each >= min_neighbours_per_type neighbours.
#   "low_entropy":    entropy in the lowest low_entropy_quantile of SCLC-A cells
#                     (most neighbours are one type).
# Both: every neighbour type with >= min_neighbours_per_type neighbours must have
# >= min_reference_cells reference cells (enough to train its VAE), and no skip_labels.
# Candidates are ranked by how close their neighbour % is to the average of the group they
# come from (all SCLC-A cells for "representative", low-entropy cells for "low_entropy");
# the closest is the query cell.
#
# Output (results/vae_steps/<query cell>/):
#   1_average_profile.csv   average % of neighbours per type, over the group
#   1_candidates.csv        every candidate, closest to the average first
#   1_query_cell.csv        the chosen cell: one row per neighbour type
#
# Run: python 1_pick_query_cell.py

import json
import os

import numpy as np
import pandas as pd

pd.set_option("display.width", 200)


# --------------------------------------------------
# Settings
# --------------------------------------------------

root = "/Users/shivanim/Documents/peer_lab/scimilarity_test"
neighbours_csv = os.path.join(root, "results/annotation_entropy_vary_k/"
                                    "annotation_entropy_all_true_counts_k50.csv")
labels_tsv = os.path.join(root, "models/model_v1.1/annotation/reference_labels.tsv")
steps_dir = os.path.join(root, "results/vae_steps")

pick = "low_entropy"         # "representative" or "low_entropy"
query_group = "SCLC-A"
k = 50
n_types = 5                  # "representative": top types that must pass the checks
min_reference_cells = 2000   # fewer than this -> too few cells to train a VAE
min_neighbours_per_type = 2  # types with fewer neighbours (1 of 50) are ignored
entropy_window = 0.15        # "representative": entropy within SCLC-A median +- this
low_entropy_quantile = 0.10  # "low_entropy": lowest 10% of SCLC-A entropy
skip_labels = ["native cell"]  # placeholder label, not a real cell type


# --------------------------------------------------
# Step 1. Neighbours of every SCLC-A cell
# --------------------------------------------------

nbrs = pd.read_csv(neighbours_csv, index_col=0)
nbrs = nbrs[nbrs["cell_type"] == query_group]
# rows = cells, columns = reference types, values = % of the 50 neighbours
pct = pd.DataFrame([json.loads(h) for h in nbrs["hits"]], index=nbrs.index).fillna(0) * 100 / k
entropy = nbrs["annotation_entropy"]

ref_counts = pd.read_csv(labels_tsv, sep="\t", header=None, usecols=[0])[0].value_counts()
print(f"{len(pct)} {query_group} cells, median annotation entropy {entropy.median():.2f}")


# --------------------------------------------------
# Step 2. Group to pick from, and its average neighbour profile
# --------------------------------------------------

if pick == "representative":
    lo, hi = entropy.median() - entropy_window, entropy.median() + entropy_window
    group = pct.index[(entropy >= lo) & (entropy <= hi)]
    avg = pct.mean()  # average over all SCLC-A cells
    print(f"Picking a representative cell: entropy {lo:.2f}-{hi:.2f}")
elif pick == "low_entropy":
    hi = entropy.quantile(low_entropy_quantile)
    group = pct.index[entropy <= hi]
    avg = pct.loc[group].mean()  # average over the low-entropy cells
    print(f"Picking a low-entropy cell: lowest {100 * low_entropy_quantile:.0f}% of "
          f"{query_group} cells, entropy <= {hi:.2f} ({len(group)} cells)")
else:
    raise ValueError(f"pick must be 'representative' or 'low_entropy', not {pick!r}")

avg = avg.sort_values(ascending=False)
avg_table = pd.DataFrame({"avg_pct_of_neighbours": avg.round(2),
                          "reference_cells": ref_counts[avg.index]})
avg_table.index.name = "type"
print(f"\nAverage neighbour profile of that group (top 10 types):")
print(avg_table.head(10).to_string())


# --------------------------------------------------
# Step 3. Candidates, ranked by closeness to the average profile
# --------------------------------------------------

min_pct = 100 * min_neighbours_per_type / k
rows = []
for cell in group:
    p = pct.loc[cell].sort_values(ascending=False)
    vae_types = p.index[p >= min_pct]  # types that would get a VAE
    if pick == "representative" and (len(vae_types) < n_types or (p.iloc[:n_types] < min_pct).any()):
        continue
    if any(t in skip_labels for t in vae_types) or (ref_counts[vae_types] < min_reference_cells).any():
        continue
    row = {"cell": cell, "entropy": round(entropy[cell], 2),
           "distance_to_average": round(np.sqrt(((pct.loc[cell] - avg) ** 2).sum()), 1),
           "n_vae_types": len(vae_types), "vae_types_cover_pct": p[vae_types].sum()}
    for i, t in enumerate(p.index[:n_types], 1):
        if p[t] > 0:
            row[f"type{i}"] = f"{t} ({p[t]:.0f}%)"
    rows.append(row)

cand = pd.DataFrame(rows).sort_values("distance_to_average").reset_index(drop=True)
group_dist = np.sqrt(((pct.loc[group] - avg) ** 2).sum(axis=1))
cand["closer_than_pct_of_group"] = [round(100 * (group_dist > d).mean(), 1)
                                    for d in cand["distance_to_average"]]
print(f"\nCandidates: {len(cand)} cells (every neighbour type with >= {min_neighbours_per_type} "
      f"neighbours has >= {min_reference_cells} reference cells). Closest 5 to the average:")
print(cand.head(5).drop(columns=["closer_than_pct_of_group"]).to_string(index=False))


# --------------------------------------------------
# Step 4. The query cell: save in its own folder
# --------------------------------------------------

query_cell = cand.loc[0, "cell"]
cell_dir = os.path.join(steps_dir, query_cell)
os.makedirs(cell_dir, exist_ok=True)

p = pct.loc[query_cell]
p = p[p > 0].sort_values(ascending=False)
q = pd.DataFrame({"type": p.index, "pct_of_neighbours": p.values,
                  "n_neighbours": (p.values * k / 100).round().astype(int),
                  "reference_cells": ref_counts[p.index].values,
                  "train_vae": p.values >= min_pct})
q.insert(0, "cell", query_cell)
q.to_csv(os.path.join(cell_dir, "1_query_cell.csv"), index=False)
avg_table.to_csv(os.path.join(cell_dir, "1_average_profile.csv"))
cand.to_csv(os.path.join(cell_dir, "1_candidates.csv"), index=False)
with open(os.path.join(steps_dir, "current_query_cell.txt"), "w") as f:
    f.write(query_cell + "\n")

print(f"\nQuery cell ({pick}): {query_cell}")
print(f"  annotation entropy {entropy[query_cell]:.2f} ({query_group} median {entropy.median():.2f})")
print(f"  neighbour mix closer to the group average than "
      f"{cand.loc[0, 'closer_than_pct_of_group']}% of the group")
print(q.drop(columns="cell").to_string(index=False))
print("\nSaved in", cell_dir)
print("Steps 2-4 will use this cell (written to current_query_cell.txt)")
