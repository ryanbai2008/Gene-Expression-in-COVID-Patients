import os
os.environ['OMP_NUM_THREADS'] = '1'

import pandas as pd
import numpy as np
import random
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import kneighbors_graph
import scipy.sparse as sparse

import igraph as ig
import leidenalg

N_NEIGHBORS = 10
RESOLUTION_PARAMETER = 1.0
RANDOM_SEED = 42

data_df = pd.read_csv("normalized_gene_data_log2.csv", index_col=False)
data_df = data_df.T

data_df = data_df.iloc[1:, :]


# Standardize features so each gene contributes equally.
scaler = StandardScaler()
X = scaler.fit_transform(data_df)

# Build a k-nearest neighbor graph based on Euclidean distance.
knn = kneighbors_graph(X, n_neighbors=N_NEIGHBORS, mode="distance", include_self=False, metric="euclidean")

# Convert distances into similarity weights for Leiden using sparse operations.
knn.data = 1.0 / (1.0 + knn.data)
knn.eliminate_zeros()
weights = knn.maximum(knn.T).tocoo()
weights = sparse.triu(weights, k=1).tocoo()

# Create an undirected igraph graph from the sparse kNN weights.
edges = list(zip(weights.row.tolist(), weights.col.tolist()))
edge_weights = weights.data.tolist()

graph = ig.Graph(edges=edges, directed=False)
graph.es["weight"] = edge_weights

graph_density = graph.density()
degrees = np.array(graph.degree())
print("Graph diagnostics:")
print(f"  n_neighbors = {N_NEIGHBORS}")
print(f"  nodes = {graph.vcount()}")
print(f"  edges = {graph.ecount()}")
print(f"  density = {graph_density:.4f}")
print(f"  avg degree = {degrees.mean():.2f}")
print(f"  min degree = {degrees.min():.0f}")
print(f"  max degree = {degrees.max():.0f}")

if graph_density < 0.02:
    print("  Graph is sparse. Consider increasing n_neighbors to add more edges.")
elif graph_density > 0.10:
    print("  Graph is dense. Consider decreasing n_neighbors to focus on stronger local structure.")
else:
    print("  Graph density is in a moderate range for Leiden clustering.")

print("Parameter guidance:")
print("  - Increase RESOLUTION_PARAMETER for more clusters.")
print("  - Decrease RESOLUTION_PARAMETER for fewer clusters.")
print("  - Adjust N_NEIGHBORS to change graph connectivity, not final k.")
print("")

# Run Leiden algorithm.
partition = leidenalg.find_partition(
    graph,
    leidenalg.RBConfigurationVertexPartition,
    weights="weight",
    resolution_parameter=RESOLUTION_PARAMETER,
    seed=RANDOM_SEED,
)

cluster_labels = np.array(partition.membership) + 1
result_df = pd.DataFrame({"sample": data_df.index, "cluster": cluster_labels})
result_df.to_csv(os.path.join(os.path.dirname(__file__), "leiden_clusters.csv"), index=False)
print(result_df.head())

# Create a color map for each cluster so nodes are colored by cluster.
cluster_ids = np.unique(cluster_labels)
palette = ["red", "blue", "green", "orange", "purple", "cyan", "magenta", "brown", "olive", "darkblue", "gold", "pink"]
cluster_color_map = {cluster_id: palette[i % len(palette)] for i, cluster_id in enumerate(cluster_ids)}
vertex_colors = [cluster_color_map[label] for label in cluster_labels]

random.seed(RANDOM_SEED)
ig.set_random_number_generator(random)

# Plot the Leiden k-NN graph and save it to a PNG file.
layout = graph.layout_fruchterman_reingold()
plot_path = os.path.join(os.path.dirname(__file__), "leiden_graph.png")
visual_style = {
    "layout": layout,
    "vertex_size": 20,
    "vertex_color": vertex_colors,
    "vertex_label": None,
    "edge_width": 0.2,
    "bbox": (1000, 1000),
    "margin": 50,
}



ig.plot(graph, target=plot_path, **visual_style)
print(f"Saved Leiden graph image to {plot_path}")
print("Saved Leiden cluster labels to clustering/leiden_clusters.csv")
