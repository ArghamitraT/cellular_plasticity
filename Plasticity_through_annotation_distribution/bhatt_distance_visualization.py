import os
import pandas as pd
import utils_AT
import scipy.cluster.hierarchy as sch
import seaborn as sns
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go



def plot_interactive_graph(G, file_name):
    pos = nx.spring_layout(G)  # or any other layout

    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            colorscale='YlGnBu',
            size=10,
            line=dict(width=2)),
        text=[str(node) for node in G.nodes()],  # Define the node labels
        hovertext=[str(node) for node in G.nodes()]
    )

    # Create a trace for text labels separate from the node markers
    label_trace = go.Scatter(
        x=node_x,
        y=node_y,
        text=[str(node) for node in G.nodes()],
        mode='text',
        hoverinfo='none',
        textposition="middle right",
        textfont=dict(
            family="Arial",
            size=12,
            color="#000000"
        )
    )

    fig = go.Figure(data=[edge_trace, node_trace, label_trace],  # Add label_trace to the data
                    layout=go.Layout(
                        showlegend=False,
                        hovermode='closest',
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        margin=dict(b=20,l=5,r=5,t=40),
                        title=dict(text=file_name, xref="paper", yref="paper", x=0.5, y=0.95)
                    ))

    fig.write_html(file_name)
    fig.show()

# Assuming G is your graph
# plot_interactive_graph(G)



# Assuming G is your graph
# plot_interactive_graph(G)

# def plot_interactive_graph(G):
#     pos = nx.spring_layout(G)  # or any other layout
#
#     edge_x = []
#     edge_y = []
#     for edge in G.edges():
#         x0, y0 = pos[edge[0]]
#         x1, y1 = pos[edge[1]]
#         edge_x.extend([x0, x1, None])
#         edge_y.extend([y0, y1, None])
#
#     edge_trace = go.Scatter(
#         x=edge_x, y=edge_y,
#         line=dict(width=0.5, color='#888'),
#         hoverinfo='none',
#         mode='lines')
#
#     node_x = []
#     node_y = []
#     for node in G.nodes():
#         x, y = pos[node]
#         node_x.append(x)
#         node_y.append(y)
#
#     node_trace = go.Scatter(
#         x=node_x, y=node_y,
#         mode='markers',
#         hoverinfo='text',
#         marker=dict(
#             showscale=True,
#             colorscale='YlGnBu',
#             size=10,
#             line=dict(width=2)))
#
#     node_adjacencies = []
#     node_text = []
#     for node, adjacencies in enumerate(G.adjacency()):
#         node_adjacencies.append(len(adjacencies[1]))
#         node_text.append(str(node))
#
#     node_trace.marker.color = node_adjacencies
#     node_trace.hovertext = node_text
#
#     # Separate trace for text to ensure better control
#     text_trace = go.Scatter(
#         x=node_x,
#         y=node_y,
#         text=node_text,
#         mode='text',
#         hoverinfo='none',
#         textposition="top center",
#         textfont=dict(
#             family="sans serif",
#             size=12,
#             color="#000000"
#         )
#     )
#
#     fig = go.Figure(data=[edge_trace, node_trace, text_trace],  # Include text trace here
#                     layout=go.Layout(
#                         title='<br>Network graph made with Python',
#                         titlefont_size=16,
#                         showlegend=False,
#                         hovermode='closest',
#                         margin=dict(b=20, l=5, r=5, t=40),
#                         annotations=[dict(
#                             text="Python code: <a href='https://www.plotly.com/'> Plotly</a>",
#                             showarrow=False,
#                             xref="paper", yref="paper",
#                             x=0.005, y=-0.002)],
#                         xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
#                         yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
#                     )
#
#     fig.show()

# Assuming G is your graph
# plot_interactive_graph(G)



def create_full_graph(dist_matrix, labels):
    G = nx.Graph()

    # Add nodes
    for label in labels:
        G.add_node(label)

    # Add all edges
    for i, label1 in enumerate(labels):
        for j, label2 in enumerate(labels):
            if i != j:  # Exclude self-loops
                # Use 1/distance to represent edge weight (optional)
                G.add_edge(label1, label2, weight=1/dist_matrix[i][j])

    return G

def plot_graph(G, title):
    plt.figure(figsize=(12, 12))
    pos = nx.spring_layout(G, weight='weight')  # Position nodes using Fruchterman-Reingold force-directed algorithm
    edge_weights = nx.get_edge_attributes(G, 'weight')
    nx.draw(G, pos, with_labels=True, node_size=500, node_color='lightblue', font_size=10, font_weight='bold',
            edge_color='gray', width=[w*10 for w in edge_weights.values()])
    plt.title(title)
    plt.savefig("../figures/"+utils_AT.create_image_name(title))
    plt.show()


def create_graph_closest_neighbors(dist_matrix, labels, num_neighbors=3):
    G = nx.Graph()

    # Add nodes
    for label in labels:
        G.add_node(label)

    # Add edges for closest neighbors
    for i, label1 in enumerate(labels):
        # Exclude the distance to the same cell type
        distances = list(dist_matrix[i])
        distances[i] = np.inf

        # Find indices of the closest neighbors
        closest_indices = np.argsort(distances)[:num_neighbors]

        for j in closest_indices:
            label2 = labels[j]
            G.add_edge(label1, label2, weight=dist_matrix[i][j])

    return G


import matplotlib.pyplot as plt
import networkx as nx


def plot_graph_clean(G, title):
    # plt.figure(figsize=(20, 20))  # Increased figure size
    # pos = nx.spring_layout(G)  # Position nodes using Fruchterman-Reingold force-directed algorithm
    #
    # # Draw nodes
    # nx.draw_networkx_nodes(G, pos, node_size=700, node_color='lightblue', alpha=0.6)
    #
    # # Draw edges
    # nx.draw_networkx_edges(G, pos, alpha=0.2)
    #
    # # Draw node labels
    # nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold')

    plt.figure(figsize=(20, 20))  # Increased figure size
    pos = nx.spring_layout(G, k=0.1, iterations=50)  # Experiment with k (optimal distance) and iterations

    # Draw the network
    nx.draw_networkx_edges(G, pos, alpha=0.2)
    nx.draw_networkx_nodes(G, pos, node_size=700, node_color='lightblue', alpha=0.6)

    # Draw labels with adjusted positions to reduce overlap
    label_pos = {k: [v[0], v[1] + 0.02] for k, v in pos.items()}  # Adjust label position
    nx.draw_networkx_labels(G, label_pos, font_size=8, font_weight='bold', alpha=0.7)

    # Set the plot title
    plt.title(title, size=20)

    # Remove the axes
    plt.axis('off')
    plt.savefig("../figures/" + utils_AT.create_image_name(title))
    # Show plot
    plt.show()



# def plot_graph(G):
#     plt.figure(figsize=(12, 12))
#     pos = nx.spring_layout(G, weight='weight')  # Position nodes using Fruchterman-Reingold force-directed algorithm
#     edge_weights = nx.get_edge_attributes(G, 'weight')
#     nx.draw(G, pos, with_labels=True, node_size=500, node_color='lightblue', font_size=10, font_weight='bold',
#             edge_color='gray', width=[w * 10 for w in edge_weights.values()])
#     plt.title("Graph of Cell Types with 3 Closest Neighbors")
#     plt.show()


def plot_heatmap_with_dendrogram(dist_matrix, labels, title):
    # Perform hierarchical clustering
    linkage = sch.linkage(sch.distance.squareform(dist_matrix), method='average')

    # Create a dendrogram
    dendro = sch.dendrogram(linkage, labels=labels, leaf_rotation=90)

    # Reorder the distances matrix and labels according to the dendrogram
    ordered_labels = [labels[i] for i in dendro['leaves']]
    ordered_matrix = dist_matrix[:, dendro['leaves']][dendro['leaves'], :]

    # Plot the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(ordered_matrix, cmap='viridis', xticklabels=ordered_labels, yticklabels=ordered_labels)
    plt.title(title)
    plt.show()

def bhattacharyya_distance(mean1, cov1, mean2, cov2):
    # Calculate the average covariance matrix
    Sigma = 0.5 * (cov1 + cov2)

    mean1 = mean1.squeeze()
    mean2 = mean2.squeeze()

    # Compute the Bhattacharyya distance
    diff_mean = mean2 - mean1
    term1 = 0.125 * np.dot(np.dot(diff_mean.T, np.linalg.inv(Sigma)), diff_mean)
    # term2 = 0.5 * np.log(np.linalg.det(Sigma) / np.sqrt(np.linalg.det(cov1) * np.linalg.det(cov2)))
    sign1, logdet1 = np.linalg.slogdet(cov1)
    sign2, logdet2 = np.linalg.slogdet(cov2)
    sign, logdet = np.linalg.slogdet(Sigma)

    # Make sure that all determinants are positive (sign should be 1)
    if sign1 == sign2 == sign == 1:
        term2 = 0.5 * (logdet - 0.5 * (logdet1 + logdet2))
    else:
        # Handle the case where the determinant is negative or zero
        # This might indicate an issue with the covariance matrices
        print("Error: Non-positive determinant encountered.")

    print("running")
    return term1 + term2


#def kl_divergence(mu1, sigma1, mu2, sigma2):
    # [The rest of the function remains the same]

def compute_distances(folder_path, cell_types):
    # Initialize a matrix to store the distances
    num_cell_types = len(cell_types)
    distances = np.zeros((num_cell_types, num_cell_types))

    # Compute pairwise distances
    for i, type1 in enumerate(cell_types):
        mean1 = np.load(os.path.join(folder_path, f'{type1}_mean.npy')).squeeze()
        cov1 = np.load(os.path.join(folder_path, f'{type1}_cov.npy')).squeeze()

        for j, type2 in enumerate(cell_types):
            if j > i:  # To avoid redundant computations
                mean2 = np.load(os.path.join(folder_path, f'{type2}_mean.npy')).squeeze()
                cov2 = np.load(os.path.join(folder_path, f'{type2}_cov.npy')).squeeze()

                distance = bhattacharyya_distance(mean1, cov1, mean2, cov2)
                # Or use KL divergence if preferred
                # distance = kl_divergence(mean1, cov1, mean2, cov2)

                distances[i, j] = distance
                distances[j, i] = distance  # Since the distance is symmetric

    return distances

# Folder paths
# folder_path1 = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/fitted_gaussian')
# folder_path2 = os.path.join(os.getcwd(), '../../data/scimilarity_data_division/fitted_gaussian_noGMM')

cell_type_data = pd.read_csv(os.path.join(os.getcwd(), '../../files/scimilarity_trainingset_cellnum_AuthorLable_modified.csv'))
cell_types = cell_type_data['Cell_Type'].tolist()  # Replace 'Cell_Type' with the actual column name

# Compute distances for each folder
distances_WGMM = np.load("bhatt_dist_W_GMM.npy")
distances_NoGMM = np.load("bhatt_dist_NO_GMM.npy")

# Assuming distances_folder1 is your distance matrix and cell_types is your list of cell types
# plot_heatmap_with_dendrogram(distances_folder1, cell_types, 'Bhatt dist Gaussian with GMM')
# plot_heatmap_with_dendrogram(distances_folder2, cell_types, 'Bhatt dist Gaussian NO GMM')

# Assuming distances_folder1 is your distance matrix and cell_types is your list of cell types
#G = create_full_graph(distances_WGMM, cell_types)
#plot_graph(G, "closest_3_neighbors_graph_")

G = create_graph_closest_neighbors(distances_WGMM, cell_types)
plot_interactive_graph(G, "../figures/bhatt_dist_W_GMM.html")

G2 = create_graph_closest_neighbors(distances_NoGMM, cell_types)
plot_interactive_graph(G2, "../figures/bhatt_dist_No_GMM.html")
#plot_graph_clean(G, "closest_3_neighbors_graph_")
print()

#
#
# plt.figure(figsize=(12, 6))
#
# # Heatmap for the first matrix
# plt.subplot(1, 2, 1)  # 1 row, 2 columns, 1st subplot
# plt.imshow(distances_folder1, cmap='viridis')
# plt.colorbar()
# plt.title('Bhatt dist Gaussian with GMM')
#
# # Heatmap for the second matrix
# plt.subplot(1, 2, 2)  # 1 row, 2 columns, 2nd subplot
# plt.imshow(distances_folder2, cmap='viridis')
# plt.colorbar()
# plt.title('Bhatt dist Gaussian no GMM')
#
# plt.savefig(utils_AT.create_image_name("Bhatt_dist_"))
# plt.show()
#
# plt.figure(figsize=(6, 6))
# plt.imshow(distances_folder1, cmap='viridis', alpha=0.5)  # alpha for transparency
# plt.imshow(distances_folder2, cmap='magma', alpha=0.5)  # different colormap
# plt.colorbar()
# plt.title('Overlay of bhatt_dist with GMM and no GMM')
# plt.savefig(utils_AT.create_image_name("Bhatt_dist_overlay"))
# plt.show()
# print()
# If you want to do something with distances_folder1 and distances_folder2, you can do it here.
