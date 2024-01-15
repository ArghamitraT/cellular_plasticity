# I believe this file is simulating plasticity (AT)
# Probably if I could get the environment, could save some time.
import numpy as np
from scipy.stats import multivariate_normal
import pandas as pd
from tqdm.auto import tqdm
import networkx as nx

import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
sns.set_context('talk')
sns.set_style('white')

import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'


from anndata import AnnData
import scanpy as sc

# Import icecream for debugging
from icecream import ic

# %%
#%load_ext autoreload
#%autoreload 2

# %%
n_dim = 15
random_seed = 45653
dens_decay = .9
var_decay = 1.5
curvature = .2
sample_res = 20

# branching structure (cells per sub-branch):
n_samples = (5*sample_res, [
    (10*sample_res, [
        (3*sample_res, [
             (10*sample_res, []), (6*sample_res, [(20*sample_res, []),]),
        ]),
        (1*sample_res, [
            (7*sample_res, [(2*sample_res, [(10*sample_res, [(1*sample_res, [(7*sample_res, []),]),]),]),]),
        ]),
    ]),
    (6*sample_res, [
        (5*sample_res, []), (8*sample_res, [(12*sample_res, []), (5*sample_res, []), ]),
    ]),
])

# %%
#%%time
np.random.seed(random_seed)

# (AT) try to understand what is happening
def sample_branch(base, velocity, sample_struck, curvature=.5, var_decay=1, dens_decay=1, n_dim=10, branch_name='b'):
    n_draws, sub_struck = sample_struck
    scale = np.sqrt(np.sum(velocity**2))
    delta = np.random.normal(0, scale, n_dim)
    new_velo = (velocity*(1-curvature) + delta*curvature) * dens_decay
    mu = base + (new_velo/2)
    q, r = np.linalg.qr(new_velo[:, None], mode='complete')
    mr = 5e-1 * np.exp(-np.arange(n_dim)*var_decay/2) * np.abs(r[0])
    mr = np.clip(mr, 1e-4*np.abs(r[0]), None)
    bases = q*mr[None, :]
    cov = bases.dot(bases.T)
    dist = multivariate_normal(mu, cov)
    samples = dist.rvs(n_draws)
    samp_list, dist_list, draws_list = [samples, ], [dist, ], [n_draws, ]
    name_list = [np.repeat(branch_name, n_draws), ]
    new_base = base + new_velo
    for i, ss in enumerate(sub_struck):
        samples, dist, n_draws, names = sample_branch(
            new_base, new_velo, ss, curvature=curvature, var_decay=var_decay,
            n_dim=n_dim, branch_name=f'{branch_name}-{i}',
        )
        samp_list += samples
        dist_list += dist
        draws_list += n_draws
        name_list += names
    return samp_list, dist_list, draws_list, name_list


sample_list, distributions, draws_list, name_list = sample_branch(
    np.zeros(n_dim), np.ones(n_dim), n_samples, n_dim=n_dim,
    dens_decay=dens_decay, var_decay=var_decay, curvature=curvature,
)

# (AT) I believe here we are forming the mutation or actual branch??
samples = np.concatenate(sample_list, axis=0)
n_draws = np.array(draws_list)
weights = n_draws / np.sum(n_draws)
pdf = np.zeros(np.sum(n_draws))
for w, dist in zip(weights, distributions):
    pdf += w*dist.pdf(samples)

# %%
sim_ad = AnnData(
    samples,
    obs=pd.DataFrame({
        'ground_truth': pdf,
        'log_ground_truth': np.log10(pdf),
        'branch_name': np.concatenate(name_list, axis=0),
    }, index=range(np.sum(n_draws))),
)
sim_ad.obs['branch_name'] = sim_ad.obs['branch_name'].astype('category')
sc.pp.neighbors(sim_ad, n_neighbors=2)
dd = sim_ad.obsp['distances'].copy()
dd.setdiag(0)
dd.eliminate_zeros()
distances = dd.data
sim_ad.obs['distance_to_closest'] = distances
sim_ad.obs_names = 'Cell_' + sim_ad.obs_names.astype(str)

# %%
sim_ad.shape

# %%
#%%time
sc.pp.neighbors(sim_ad, n_neighbors=30)
sc.tl.umap(sim_ad)
sc.tl.leiden(sim_ad)
sim_ad

# %%
sc.pl.scatter(sim_ad, basis='umap', color=['branch_name'], size=10)

# %%
sim_ad

# %%
"""
# Generate phylogenetic tree
"""

# %%
dcs = sim_ad.to_df()
dcs.head()

# %%
# Compute pairwise distance matrix
from scipy.spatial.distance import pdist, squareform
dists = pd.DataFrame(squareform(pdist(dcs)))

dists.index = dcs.index
dists.columns = dcs.index
dists[dists==0] = np.nan
ic(dists.shape)

# %%
umap = pd.DataFrame(sim_ad.obsm['X_umap'], index=sim_ad.obs_names).merge(sim_ad.obs[['branch_name']], left_index=True, right_index=True)
umap.head()

# %%
# Select the rightmost cell as the root
root = umap[umap[0]==umap[0].min()].index[0]
root

# %%
# Plot the UMAP and highlight the root cell in red
plt.figure(figsize=(5,5))
plt.scatter(umap[0], umap[1], c='grey', s=0.5)
plt.scatter(umap.loc[root, 0], umap.loc[root, 1], c='red', s=10)
plt.axis('off')
plt.show()

# %%


# %%
# Compute tree by neighbor joining
import phylo_simlib as phylo
t = phylo.neighbor_joining(dists, outgroup=root)

# %%
t.write(outfile='small_simbranch_tree.nw')
sim_ad.write('small_simbranch.h5ad')


# %%
tree_dm = phylo.get_distance_matrix_from_tree(t)
# Ensure the leaves are in the same order as the anndata object
tree_dm = tree_dm.loc[sim_ad.obs_names, sim_ad.obs_names]

# %%
diff_dm = dists
diff_dm.fillna(0, inplace=True)

# %%
x = np.triu(diff_dm.values).flatten()
y = np.triu(tree_dm.values).flatten()

# %%
plt.figure(figsize=(5, 5))
plt.scatter(x,y, s=1, alpha=.5)
plt.xlabel('Phenotypic Distance')
plt.ylabel('Tree Distance')
# Annotate with correlation
plt.text(0.1, 0.9, f'Pearson R = {np.corrcoef(np.triu(diff_dm.values).flatten(), np.triu(tree_dm.values).flatten())[0,1]:.2f}', transform=plt.gca().transAxes)

# %%
# Make a copy of the distance matrix
leiden_pairs = diff_dm.copy()
for col in leiden_pairs.columns:
    leiden_pairs[col] = sim_ad.obs['leiden'].loc[col] + ' - ' + sim_ad.obs['leiden'].astype(str)
leiden_pairs.head()

# %%
# Replace pairs in {leiden1}-{leiden2} format so that leiden1 < leiden2
for col in leiden_pairs.columns:
    leiden_pairs[col] = leiden_pairs[col].apply(lambda x: '-'.join([str(i) for i in sorted([int(i) for i in x.split('-')])]))
leiden_pairs.head()

# %%
unique_leiden_pairs = np.unique(leiden_pairs.values.flatten())

# Generate a color map for the leiden pairs
leiden_cmap = dict(zip(unique_leiden_pairs, sns.color_palette('tab20', len(unique_leiden_pairs))))
leiden_annots = np.triu(leiden_pairs.values).flatten()

# %%
leiden_cmap

# %%
# Plot x vs y, colored by leiden pairs
plt.figure(figsize=(5, 5))
for pair in unique_leiden_pairs:
    plt.scatter(x[leiden_annots==pair], y[leiden_annots==pair], s=1, alpha=.5, color=leiden_cmap[pair], label=pair)
plt.xlabel('Phenotypic Distance')
plt.ylabel('Tree Distance')
# Plot legend in multiple columns
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', ncol=5)
# Annotate with correlation
plt.text(0.1, 0.9, f'Pearson R = {np.corrcoef(np.triu(diff_dm.values).flatten(), np.triu(tree_dm.values).flatten())[0,1]:.2f}', transform=plt.gca().transAxes)

# %%
# Plot x vs y for each leiden pair
for pair in unique_leiden_pairs:
    if np.corrcoef(x[leiden_annots==pair], y[leiden_annots==pair])[0,1] > 0.7:
        continue

    plt.figure(figsize=(2,2))
    plt.scatter(x[leiden_annots==pair], y[leiden_annots==pair], s=1, alpha=.5, color=leiden_cmap[pair])
    plt.xlabel('Phenotypic Distance')
    plt.ylabel('Tree Distance')
    plt.title(f'Correlation for pair {pair}: {np.corrcoef(x[leiden_annots==pair], y[leiden_annots==pair])[0,1]:.2f}')
    plt.show()

# %%


# %%
# Plot anndata colored by leiden with labels on UMAP
with plt.rc_context({'figure.figsize': (10, 10)}):
    sc.pl.umap(sim_ad, color=['leiden'], size=20, legend_loc='on data')

# %%
# Choose a random cell from leiden cluster 8 and plot the umap coloured by phenotypic distance to that cell
leiden_8 = np.random.choice(sim_ad.obs_names[sim_ad.obs['leiden']=='9'], size=1)
sim_ad.obs['pheno_dist/tree_dist'] = diff_dm[leiden_8] / tree_dm[leiden_8]
with plt.rc_context({'figure.figsize': (10, 10)}):
    sc.pl.umap(sim_ad, color=['pheno_dist/tree_dist'], size=20)


# %%
diff_dm[leiden_8] / tree_dm[leiden_8]

# %%
leiden_8

# %%
"""
# Overlay CRISPR mutations onto the tree
We use Cassiopeia simulation engine to overlay CRISPR indels on the tree.

"""

# %%
# We need to convert the tree to a CassiopeiaTree object first

# %%
from cassiopeia.data import CassiopeiaTree
cass_tree = CassiopeiaTree()
cass_tree.populate_tree(t)

assert set(cass_tree.leaves) == set(sim_ad.obs_names)

# %%
# instantiate Cas9 lineage tracing object & overlay data onto ground_truth_tree
# We will use the default parameters for the simulator
# See documentation for more details: https://cassiopeia-lineage.readthedocs.io/en/latest/api/reference/cassiopeia.sim.Cas9LineageTracingDataSimulator.html#cassiopeia.sim.Cas9LineageTracingDataSimulator

from cassiopeia.sim import Cas9LineageTracingDataSimulator

lt_sim = Cas9LineageTracingDataSimulator(number_of_cassettes=10,
                                        size_of_cassette=3,
                                        mutation_rate=0.01,
                                        number_of_states=100,
                                        state_priors=None,
                                        heritable_silencing_rate=0.0001,
                                        stochastic_silencing_rate=0.01,
                                        heritable_missing_data_state=-1,
                                        stochastic_missing_data_state=-1,
                                        random_seed=None,
                                        collapse_sites_on_cassette=True)
lt_sim.overlay_data(cass_tree)

# %%
cass_tree.character_matrix

# %%
sim_ad.shape

# %%
"""
# Save simulated data to output files
"""

# %%
# We need to (1) pickle the CassiopeiaTree object, (2) save the anndata object, (3) save the distance matrices, (4) save the tree in newick format and (5) save the character matrix

outdir = 'for_yang_evaluation_plasticity/'
import os
os.makedirs(outdir, exist_ok=True)

# (1) Pickle the CassiopeiaTree object
import pickle
with open(os.path.join(outdir, 'cass_tree.pkl'), 'wb') as f:
    pickle.dump(cass_tree, f)

# (2) Save the anndata object
sim_ad.write(os.path.join(outdir, 'sim_ad.h5ad'))

# (3) Save the distance matrices
diff_dm.to_csv(os.path.join(outdir, 'diff_dm.csv'))
tree_dm.to_csv(os.path.join(outdir, 'tree_dm.csv'))

# (4) Save the tree in newick format
t.write(outfile=os.path.join(outdir, 'tree.nw'))

# (5) Save the character matrix
cass_tree.character_matrix.to_csv(os.path.join(outdir, 'character_matrix.csv'))

# %%
