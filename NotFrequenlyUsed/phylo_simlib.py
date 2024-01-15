import pandas as pd
import numpy as np
from icecream import ic
import networkx as nx


def neighbor_joining(distance_matrix, outgroup=None):
    """
    Neighbor joining algorithm for constructing a phylogenetic tree from a distance matrix.
    :param: distance_matrix - pd.DataFrame containing inter-leaf distances
    :param: outgroup - name of the outgroup leaf
    :return: ete3.Tree object
    """

    from ete3 import TreeNode

    leaves = distance_matrix.index
    n = len(leaves)
    D = distance_matrix.copy()
    D.replace(0, np.nan, inplace=True)
    D = D.values

    names = list(distance_matrix.index)
    if outgroup is not None:
        assert outgroup in names, f'Outgroup {outgroup} is not a valid leaf in the distance matrix'
    joins = 0

    nodes = {}
    for name in names:
        nodes[name] = TreeNode(name=name)

    ic(D.shape[0])
    while D.shape[0] > 2:
        if D.shape[0] % 50 == 0:
            ic(D.shape[0])

        # Normalized sum of distances from node to all other points
        R = np.nansum(D, axis=0) / (D.shape[0] - 2)
        R = R.reshape(-1, 1)

        R_tile = np.tile(R, D.shape[0])

        # Compute the adjusted distance matrix
        Q = D - (R_tile + R_tile.T)

        # Select the pair of nodes with minimum distance
        i, j = np.unravel_index(np.nanargmin(Q), Q.shape)

        # Replace (i,j) with a new node k and compute branch lengths (i,k) and (j,k)
        branch_i = (D[i, j] + R[i] - R[j]) / 2
        branch_j = (D[i, j] + R[j] - R[i]) / 2
        # Compute the distance between all points to k
        D_k = np.delete((np.nansum(D[[i, j]], axis=0) - D[i, j]) / 2, [i, j])

        # Construct new distance matrix with one fewer node
        _D = np.delete(np.delete(D, [i, j], axis=1), [i, j], axis=0)
        new_D = np.zeros((_D.shape[0] + 1, _D.shape[1] + 1))
        new_D[:-1, :-1] = _D
        new_D[-1, :-1] = D_k
        new_D[:-1, -1] = D_k
        new_D[-1, -1] = np.nan

        D = new_D

        # Construct tree from joins using newick string structure
        new_name = f'({names[i]}:{branch_i[0]}, {names[j]}:{branch_j[0]})'

        names.pop(max(i, j))
        names.pop(min(i, j))
        names.append(new_name)

    # Connect the final two nodes to create full unrooted tree
    newick = f'({names[0]}:{D[0, 1] / 2}, {names[1]}:{D[0, 1] / 2});'
    from ete3 import Tree
    t = Tree(newick)

    if outgroup is None:
        # Root tree at midpoint
        root_point = t.get_midpoint_outgroup()
    else:
        root_point = outgroup

    t.set_outgroup(root_point)

    return t


def name_internal_nodes(tree):
    # Name internal nodes of a tree
    i = 1
    for node in tree.traverse():
        if not node.is_leaf():
            node.name = str(i) + "_"
            i += 1
    return None


def ete3_to_nx(tree):
    """
    Convert an ete3 Tree to a networkx graph, preserving branch lengths
    """
    name_internal_nodes(tree)

    G = nx.Graph()
    for node in tree.traverse():
        if node.is_root():
            continue
        G.add_edge(node.up.name, node.name, weight=node.dist)
    return G


def get_distance_matrix_from_tree(tree):
    """
    Compute the distance matrix from a tree. Converts the tree to a networkx graph
    and uses scipy's shortest path function to compute the distance matrix.
    :param tree: ete3 tree
    :return: distance matrix
    """
    G = ete3_to_nx(tree)
    # Compute all-pairs distances from networkx graph
    from scipy.sparse.csgraph import shortest_path
    _dm = pd.DataFrame(shortest_path(nx.to_scipy_sparse_array(G), directed=False))
    _dm.index = G.nodes
    _dm.columns = G.nodes

    # Get leaves from the tree
    leaves = [l.name for l in tree.get_leaves()]

    return _dm.loc[leaves, leaves]


def plot_tree_depth_distribution(tree):
    """
    Plots the distribution of the depth of each node in the tree
    :param tree: (ete3.Tree) an ete3 tree object
    """
    depths = [t.children[0].get_distance(x.name) for x in t.children[0].get_leaves()]
    plt.figure(figsize=(4, 4))
    sns.distplot(depths)
    plt.xlabel('Depth')
    plt.ylabel('Density')
    plt.title('Distribution of Node Depths')
    plt.show()
    plt.close()


def annotate_depths(tree):
    """
    Annotates the depth of each node in the tree. The root depth is given as 0.
    :param tree: (ete3.Tree) an ete3 tree object
    """
    for node in t.traverse("preorder"):
        if node.is_root():
            node.add_feature('depth', 0)
        else:
            node.add_feature('depth', node.up.depth + 1)
    return


def plot_generation_distribution(tree):
    """
    Plots the distribution of the generation of each leaf in the tree
    """
    try:
        leaf_depths = [n.depth for n in t.get_leaves()]
    except:
        annotate_depths()
        leaf_depths = [n.depth for n in t.get_leaves()]
    plt.figure(figsize=(4, 4))
    sns.distplot(leaf_depths)
    plt.xlabel('Generation')
    plt.ylabel('Density')
    plt.title('Distribution of Leaf Generations')
    plt.show()
    plt.close()


def plot_branch_length_distribution(tree):
    """
    Plots the distribution of the branch lengths in the tree
    :param tree: (ete3.Tree) an ete3 tree object
    """
    branch_lengths = [x.dist for x in t.children[0].get_leaves()]
    plt.figure(figsize=(4, 4))
    sns.distplot(branch_lengths)
    plt.xlabel('Branch Length')
    plt.ylabel('Density')
    plt.title('Distribution of Branch Lengths')
    plt.show()
    plt.close()
