import pandas as pd
import numpy as np
import walker
from icecream import ic

import networkx as nx
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm.notebook import tqdm


def graph_from_connectivities(adj_matrix, cell_names):
    """
    Construct a networkx graph from the given adjacency matrix. Label nodes
    according to cell_names; assumes the order of nodes in the adjacency matrix
    matches the order of nodes in cell_names.

    @param: adj_matrix - scipy.sparse.csr_matrix containing a 1 where an edge
                         connects two nodes, 0 otherwise.
    @param: cell_names - list of names of cells. Must be of the same length as adj_matrix
    """
    from scipy.sparse import csr_matrix
    if isinstance(adj_matrix, csr_matrix):
        H = nx.from_scipy_sparse_array(adj_matrix)
        nx.relabel_nodes(H, dict(zip(H.nodes(), cell_names)), copy=False)
        return H

    else:
        raise AttributeError(f'''graph_from_connectivities is not implemented 
        for data type {type(adj_matrix)}''')


def add_random_edges(G, add_edge_prob):
    """
    Introduce random connections between edges in G. If the total number of potential edges
    in a fully connected version of G is |E| then add_edge_prob \times |E| edges are introduced.
    Returns a new graph containing all the original edges in G as well as the new randomly
    drawn edges.
    :param G: networkx graph
    :param add_edge_prob: float between (0,1) representing the proportion of total potential edges added.
    :return: networkx graph containing all the original edges in G as well as the new randomly drawn edges.
    """

    nodes = list(G.nodes)

    random_G = np.random.choice([0, 1],
                                size=(len(nodes), len(nodes)),
                                p=[1 - add_edge_prob / 2, add_edge_prob / 2])

    num_new_edges = sum(random_G != 0)
    ic(num_new_edges)

    # Symmetrize the matrix
    random_G = ((random_G + random_G.T) > 0).astype(int)

    random_G = pd.DataFrame(random_G)
    random_G.columns = nodes
    random_G.index = nodes

    random_G = nx.from_numpy_array(random_G)
    nx.relabel_nodes(random_G, dict(zip(random_G.nodes(), nodes)), copy=False)

    return nx.compose(G, random_G)


def get_distances_of_moves(G, sources, targets):
    """
    Computes the shortest graph distance from each source in sources
    and the corresponding target in targets. Return these disances in
    a list.
    :param G: networkx graph
    :param sources: list of source nodes
    :param targets: list of target nodes
    :return: list of distances between each source and target
    """

    assert len(sources) == len(targets), """`sources` and `targets` must be lists
                                        of the same length"""

    dists = []
    for source, target in list(zip(sources, targets)):
        dists.append(len(nx.shortest_path(G, source=source, target=target)))
    return np.array(dists)


def replace_df_rows(df, source_cells, target_cells):
    """
    Replace dataframe rows of source cells with corresponding target cells. Entries
    in source_cells and in target_cells must all be contained in the index of
    df.
    :param: df - pd.DataFrame
    :param: source_cells - list of cell names of original entries
    :param: target_cells - list of cell names which original entries should be
                           replaced with.
    """

    assert len(source_cells) == len(target_cells), '''Length of `source_cells` 
    must match the length of `target_cells`'''

    df2 = df.copy()
    df2.loc[source_cells] = df.loc[target_cells].values

    return df2

def add_plasticity_to_anndata(ad, source_cells, target_cells):
    """
    Add plasticity to an anndata object by replacing the source cells with the target cells.
    ad.X and ad.obs are modified, as well as the UMAP coordinates so that the source cells now
    occupy the same position/ have the same data as the target cells.
    :param ad: anndata object
    :param source_cells: list of source cells
    :param target_cells: list of target cells
    :return: anndata object with plasticity added
    """
    sim_ad_post_switch = ad.copy()
    X = sim_ad_post_switch.to_df()
    X.loc[source_cells, :] = X.loc[target_cells, :].values
    sim_ad_post_switch.X = X.values
    del X
    sim_ad_post_switch.obs.loc[source_cells, :] = ad.obs.loc[target_cells, :].values
    sim_ad_post_switch.obs['simulated_cell_id'] = ad.obs_names
    sim_ad_post_switch.obs['is_plastic'] = 0
    sim_ad_post_switch.obs.loc[source_cells, 'is_plastic'] = 1
    sim_ad_post_switch.obs.loc[source_cells, 'simulated_cell_id'] = target_cells
    sim_ad_post_switch.obs = sim_ad_post_switch.obs.join(ad.obs, rsuffix='_pre_plasticity')

    umap = pd.DataFrame(sim_ad_post_switch.obsm['X_umap'])
    umap.index = sim_ad_post_switch.obs_names
    umap.loc[source_cells, :] = umap.loc[target_cells, :].values
    sim_ad_post_switch.obsm['X_umap'] = umap.values

    return sim_ad_post_switch

def modified_hamming_distance(x,y):
    """
    This function calculates the modified Hamming distance between two arrays;
    the distance is the average of the distances between each element in the array.
    If two elements are the same, the distance is 0; if two elements are different and
    both are non-zero, the distance is 1; if one element is zero and the other is non-zero,
    the distance is 2. We ignore sites where either element is missing (i.e. np.nan).
    :param x: (np.array)
    :param y: (np.array)
    :return:
    """
    dists = 2 * (x != y).astype(float) - np.logical_xor(x == 0, y == 0).astype(float)
    nan_mask = np.logical_or(np.isnan(x), np.isnan(y))
    dists[nan_mask] = np.nan
    return np.nanmean(dists)

def sorted_neighbors(dm):
    """
    Returns a sorted list of neighbors for each cell in the distance matrix.
    Each row in the returned matrix is a list of neighbors in order of increasing distance
    for the corresponding cell in the distance matrix.
    :param dm: pd.DataFrame, distance matrix where rows and columns are cells
    :return: pd.DataFrame, sorted list of neighbors for each cell in the distance matrix (each row is a list)
    """
    sorted_distances = dm.apply(lambda row: row.sort_values(ascending=True).index.tolist(), axis=1)
    nhoods = np.vstack(sorted_distances.values)
    nhoods = pd.DataFrame(nhoods, index=sorted_distances.index)
    return nhoods

def get_overlaps(tree_dm, pheno_dm, radii=None):
    """
    Gets the distribution of the number of cells in the intersection of the phylogenetic and phenotypic neighborhoods
    for each radius
    :param tree_dm: (pd.DataFrame) a pandas dataframe of pairwise distances between cells in phylogenetic space
    :param pheno_dm: (pd.DataFrame) a pandas dataframe of pairwise distances between cells in phenotypic space
    :param radii: (list) a list of radii to consider
    """
    if radii is None:
        radii = np.arange(1, tree_dm.shape[0]-1, 5)

    tree_dm = tree_dm.loc[:, tree_dm.index]
    pheno_dm = pheno_dm.loc[tree_dm.index, tree_dm.columns]

    tree_nhoods = sorted_neighbors(tree_dm)
    pheno_nhoods = sorted_neighbors(pheno_dm)

    overlaps = {}
    for r in tqdm(radii):
        tree_neigh = tree_nhoods.iloc[:, :r]
        phen_neigh = pheno_nhoods.iloc[:, :r]

        intersections = []
        for row1, row2 in zip(tree_neigh.values, phen_neigh.values):
            intersection = len(np.intersect1d(row1, row2))
            intersections.append(intersection)

        overlaps[r] = intersections

    overlaps = pd.DataFrame(overlaps)
    overlaps.index = tree_dm.index
    return overlaps

def gini_index(arr):
    """Compute the Gini index for a 1D array
    :param arr: A 1D array which contains values between 0 and 1
    :return: The Gini index for the array
    """
    # Ensure that the array values lie between 0 and 1
    assert np.all(arr >= 0) and np.all(arr <= 1)

    # Sort the array
    arr = np.sort(arr)
    # Get the cumulative sum
    cumsum = np.cumsum(arr)
    # Get the cumulative proportion
    cumprop = cumsum / cumsum[-1]
    # Get the Gini index
    gini = 1 - np.sum((cumprop[1:] + cumprop[:-1]) * (arr[1:] - arr[:-1]))
    return gini

def get_gini_index(overlaps):
    """Compute the Gini index for each row in a 2D array
    :param overlaps: pd.DataFrame of overlaps for each cell for each radius of interest
    :return: pd.DataFrame of the Gini index for each cell
    """
    # Normalize the overlaps so that they are proportions of the total number of cells
    overlaps = overlaps / overlaps.shape[0]
    cell_names = overlaps.index
    overlaps = overlaps.values

    # Ensure that overlaps is a 2D array
    assert len(overlaps.shape) == 2

    # Compute the Gini index for each row
    gini = np.array([gini_index(overlaps[i]) for i in range(overlaps.shape[0])])

    return pd.DataFrame(gini, index=cell_names, columns=['Gini Index'])



def generate_plasticity(ad,
                        num_neighbors,
                        walk_length=100,
                        plot_umap=False,
                        plot_distributions=False,
                        n_highways=0):
    """
    Introduce plasticity to the system by performing a random walk on the diffusion map graph.
    Original cell identities are swapped with the cell identity of the cell at the end of the random walk.

    We can tune the level of plasticity by changing the number of neighbors used to construct the diffusion map graph,
    the length of the random walk, and the number of highways added between celltypes.

    :param ad: (anndata.AnnData) an anndata object containing the diffusion map graph
    :param num_neighbors: (int) the number of neighbors to use when constructing the diffusion map graph
    :param walk_length: (int) the length of the random walk
    :param plot_umap: (bool) if True, plots the UMAP of the cells colored by celltype
    :param plot_distributions: (bool) if True, plots the distribution of the number of connections in the neighbors graph
    :param n_highways: (int) the number of highways to add between celltypes

    :return: diffusion components (pd.DataFrame) the diffusion components of the cells after plasticity is introduced
    :return: dc_distance_moved (np.ndarray) the distance moved by each cell in diffusion component space
    """
    import plasticity_simlib as plastic

    import copy
    ad = copy.deepcopy(ad)

    # Compute the diffusion map
    sc.pp.neighbors(ad, n_neighbors=num_neighbors, use_rep='DM_EigenVectors')
    if plot_umap:
        with plt.rc_context({"figure.figsize": (4, 4), "figure.dpi": (100)}):
            sc.pl.umap(ad, color='celltype', edges=True, s=10, edges_width=0.01)

        # Plot number of connections in anndata neighbors graph
        nn_graph = ad.obsp['connectivities']
        sns.displot(nn_graph.sum(axis=1))
        plt.title('Number of connections in neighbors graph')
        plt.show()

    if n_highways > 0:
        # Add highways between celltypes
        print(f'Adding {n_highways} highways between celltypes')
        # Convert the connectivities matrix to a sparse lil_matrix
        ad.uns['neighbors']['connectivities'] = ad.uns['neighbors']['connectivities'].tolil()
        import itertools
        celltype_tuples = list(itertools.combinations(['CD8 T cells Mem', 'CD8 T cells Naive', 'T cells 2'], 2))
        n = int(n_highways / len(celltype_tuples))

        for ct1, ct2 in celltype_tuples:
            # Select n cells from each celltype
            ct1_cells = ad.obs[ad.obs['celltype'].isin([ct1])].sample(n=n, replace=True).index
            ct2_cells = ad.obs[ad.obs['celltype'].isin([ct2])].sample(n=n, replace=True).index

            # Get the numerical index of ct1 cells in the original AnnData object
            ix1 = ad.obs_names.get_indexer(ct1_cells)
            # Get the numerical index of ct2 cells in the original AnnData object
            ix2 = ad.obs_names.get_indexer(ct2_cells)
            # Construct tuples from two lists of indices, ix1 and ix2
            tuples = list(itertools.product(ix1, ix2))
            # Add a connection between each tuple of cells
            for c1, c2 in tuples:
                ad.uns['neighbors']['connectivities'][c1, c2] = 1
                ad.uns['neighbors']['connectivities'][c2, c1] = 1
        # Convert the connectivities matrix back to a sparse csr_matrix
        ad.uns['neighbors']['connectivities'] = ad.uns['neighbors']['connectivities'].tocsr()

    F = plastic.graph_from_connectivities(ad.obsp['connectivities'], ad.obs_names)
    # p is the probability of returning to the previous node,
    # q controls the probability of moving to an unexplored node; smaller q means more exploration
    walks = walker.random_walks(F, n_walks=1, walk_len=walk_length, p=0.01, q=.001)
    final_nodes = np.array(F.nodes)[walks[:, -1]]
    start_nodes = np.array(F.nodes)[walks[:, 0]]

    # # Randomly choose 50% of the indices
    # idx = np.random.choice(np.arange(len(start_nodes)), size=int(len(start_nodes)/2), replace=False)
    # final_nodes[idx] = start_nodes[idx]

    # The 'final nodes' are the nodes that the walker ended up at after the walk
    # representing transitions in cell states - this corresponds to plasticity

    # Get diffusion components for each cell
    diffusion_comps = pd.DataFrame(ad.obsm['DM_EigenVectors'])
    diffusion_comps.index = ad.obs_names

    # Compute the diffusion distance between each start and final node
    dc_distance_moved = np.linalg.norm(
        diffusion_comps.loc[start_nodes].values - diffusion_comps.loc[final_nodes].values, axis=1)
    dc_distance_moved = pd.DataFrame(dc_distance_moved)
    dc_distance_moved.index = start_nodes

    # Compute the number of steps between each start and final node
    path = dict(nx.all_pairs_shortest_path(F))

    if plot_distributions:
        # Plot the distribution of diffusion distances as a violin plot
        plt.figure(figsize=(4, 4))
        sns.violinplot(dc_distance_moved)
        plt.xlabel('Diffusion Distance')
        plt.ylabel('Density')
        plt.title('Diffusion Distance Distribution between Start and Final Nodes')
        plt.show()

    # Get the diffusion components for the final nodes
    final_dcs = diffusion_comps.loc[final_nodes]

    # Replace the diffusion components for each start node with the corresponding final node
    # diffusion components
    diffusion_comps.loc[start_nodes] = final_dcs.values

    return diffusion_comps, dc_distance_moved