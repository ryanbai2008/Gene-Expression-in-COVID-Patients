import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

os.environ["OMP_NUM_THREADS"] = "2"

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import kneighbors_graph
import igraph as ig
import leidenalg
from scipy.cluster import hierarchy
import scipy.sparse as sparse

# ------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------
N_PCS = 10            # PCs used for clustering (variance-justified choice)
K = 2                 # number of clusters for KMeans / Hierarchical / Leiden target
LEIDEN_N_NEIGHBORS = 10   # best-behaved graph: moderate density, clean resolution steps
LEIDEN_RESOLUTION = 0.20  # confirmed to give exactly 2 clusters at n_neighbors=15
N_BOOT = 200          # number of bootstrap resamples
RANDOM_STATE = 42

rng = np.random.default_rng(RANDOM_STATE)

# ------------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------------
data_df = pd.read_csv("normalized_gene_data_log2.csv", index_col=0)
data_df = data_df.T  # samples x genes

print(f"Data shape (samples x genes): {data_df.shape}")

X_scaled = StandardScaler().fit_transform(data_df)

n_pcs = min(N_PCS, X_scaled.shape[0] - 1, X_scaled.shape[1])
X_pca = PCA(n_components=N_PCS, random_state=RANDOM_STATE).fit_transform(X_scaled)

n_samples = X_pca.shape[0]
print(f"PCA space shape: {X_pca.shape}")


# ------------------------------------------------------------------
# CLUSTERING FUNCTIONS (each takes an (n, n_pcs) array, returns labels)
# ------------------------------------------------------------------
def run_kmeans(X, k=K, random_state=RANDOM_STATE):
    return KMeans(n_clusters=k, random_state=random_state, n_init='auto').fit_predict(X)


def run_hierarchical(X, k=K):
    linkage_matrix = hierarchy.linkage(X, method="ward")
    labels = hierarchy.fcluster(linkage_matrix, t=k, criterion="maxclust")
    return labels

def run_leiden(X, resolution=LEIDEN_RESOLUTION, n_neighbors=LEIDEN_N_NEIGHBORS, random_state=RANDOM_STATE):
    n = X.shape[0]
 
    knn = kneighbors_graph(X, n_neighbors=n_neighbors, mode="distance", include_self=False, metric="euclidean")
    knn.data = 1.0 / (1.0 + knn.data)
    knn.eliminate_zeros()
    weights = knn.maximum(knn.T).tocoo()
    weights = sparse.triu(weights, k=1).tocoo()
 
    edges = list(zip(weights.row.tolist(), weights.col.tolist()))
    edge_weights = weights.data.tolist()
 
    graph = ig.Graph(n=n, edges=edges, directed=False)
    graph.es["weight"] = edge_weights
 
    partition = leidenalg.find_partition(
        graph, leidenalg.RBConfigurationVertexPartition,
        weights="weight", resolution_parameter=resolution, seed=random_state
    )
    return np.array(partition.membership)


METHODS = {
    'KMeans': run_kmeans,
    'Hierarchical': run_hierarchical,
    'Leiden': run_leiden,
}

# ------------------------------------------------------------------
# ORIGINAL (FULL-DATA) CLUSTERING PER METHOD
# ------------------------------------------------------------------
original_labels = {}
for name, fn in METHODS.items():
    labels = fn(X_pca)
    original_labels[name] = labels
    n_clusters_found = len(np.unique(labels))
    print(f"{name}: {n_clusters_found} clusters on full data")
    if name == 'Leiden' and n_clusters_found != K:
        print(f"  WARNING: Leiden gave {n_clusters_found} clusters, not {K}, at "
              f"resolution={LEIDEN_RESOLUTION}, n_neighbors={LEIDEN_N_NEIGHBORS}. "
              f"Re-check the resolution scan before trusting the subsampling comparison.")

# ------------------------------------------------------------------
# Subsampling STABILITY 
# ------------------------------------------------------------------

ari_scores = {name: [] for name in METHODS}
failed_runs = {name: 0 for name in METHODS}

for b in range(N_BOOT):
    boot_idx = rng.choice(n_samples, size=int(0.8 * n_samples), replace=False)    
    X_boot = X_pca[boot_idx]

    for name, fn in METHODS.items():
        try:
            boot_labels = fn(X_boot)
            # need at least 2 clusters on both sides for ARI to be meaningful
            if len(np.unique(boot_labels)) < 2:
                failed_runs[name] += 1
                continue
            orig_subset = original_labels[name][boot_idx]
            ari = adjusted_rand_score(orig_subset, boot_labels)
            ari_scores[name].append(ari)
        except Exception:
            failed_runs[name] += 1

# ------------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------------
print("\n===== Subsampling STABILITY SUMMARY =====")
summary = {}
for name in METHODS:
    scores = ari_scores[name]
    if len(scores) == 0:
        print(f"{name}: no valid subsampling runs (all failed / degenerate)")
        continue
    mean_ari = np.mean(scores)
    std_ari = np.std(scores)
    summary[name] = mean_ari
    print(f"{name}: mean ARI = {mean_ari:.4f} (std = {std_ari:.4f}), "
          f"valid runs = {len(scores)}/{N_BOOT}, failed = {failed_runs[name]}")

if summary:
    best_method = max(summary, key=summary.get)
    print(f"\nMost stable method: {best_method} (mean ARI = {summary[best_method]:.4f})")
else:
    best_method = None
    print("\nNo method produced valid subsampling results — check LEIDEN_RESOLUTION / data size.")

# ------------------------------------------------------------------
# PLOT: distribution of ARI scores per method
# ------------------------------------------------------------------
plt.figure(figsize=(8, 5))
plot_data = [ari_scores[name] for name in METHODS if len(ari_scores[name]) > 0]
plot_labels = [name for name in METHODS if len(ari_scores[name]) > 0]

if plot_data:
    plt.boxplot(plot_data, labels=plot_labels)
    plt.ylabel('Adjusted Rand Index (vs. original clustering)')
    plt.title(f'Subsampling Cluster Stability (n={N_BOOT} resamples)')
    plt.grid(True, axis='y')
    plt.tight_layout()
    plt.savefig("subsampling_stability_comparison.png", dpi=150)
    plt.show()

# ------------------------------------------------------------------
# FINAL: cluster assignments from the most stable method
# ------------------------------------------------------------------
if best_method is not None:
    final_labels = original_labels[best_method]
    out_df = pd.DataFrame({
        'sample': data_df.index,
        'cluster': final_labels
    })
    out_df.to_csv("final_cluster_assignments.csv", index=False)
    print(f"\nSaved final cluster assignments (method={best_method}) to final_cluster_assignments.csv")