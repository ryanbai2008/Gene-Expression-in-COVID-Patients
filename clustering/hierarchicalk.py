
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import pdist
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
from scipy.cluster import hierarchy
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


# X is your data matrix: shape (n_samples, n_features)
# Load your normalized gene expression data.
data_df = pd.read_csv("normalized_gene_data_log2.csv", index_col=0)
data_df = data_df.T

# data_df = data_df.iloc[1:, :]

# Standardize before clustering.
scaler = StandardScaler()
X = scaler.fit_transform(data_df)


pca = PCA(n_components=10, random_state=42)
X = pca.fit_transform(X)

explained_variance = pca.explained_variance_ratio_
cumulative_variance = np.cumsum(explained_variance)



# 1) Build the hierarchical tree using Euclidean distance
Z = linkage(X, method='ward')  # Ward uses Euclidean geometry and is common for hierarchical clustering

# 2) Try different numbers of clusters
ks = range(2, min(11, X.shape[0]))
scores = []

for k in ks:
    labels = fcluster(Z, t=k, criterion='maxclust')
    if len(np.unique(labels))>1:
        score = silhouette_score(X, labels, metric='euclidean')
        scores.append(score)
    else:
        scores.append(np.nan)


print("Silhouette scores for each k:")
for k, score in zip(ks, scores):
    print(f"  k={k}: {score:.4f}")

best_k = ks[np.argmax(scores)]
print("Best k:", best_k)
print("Best silhouette:", max(scores))

# Get cluster labels for the best k
best_labels = fcluster(Z, t=best_k, criterion='maxclust')
cluster_df = pd.DataFrame({"sample": data_df.index, "cluster": best_labels})
cluster_df.to_csv("hierarchical_clusters_bestk.csv", index=False)
print(f"Saved hierarchical cluster assignments to hierarchical_clusters_bestk.csv")

# 3) Plot silhouette vs k
plt.figure(figsize=(8, 5))
plt.plot(list(ks), scores, marker='o')
plt.xlabel("Number of clusters")
plt.ylabel("Silhouette score")
plt.title("Hierarchical clustering: choose k")
plt.grid(True)
plt.tight_layout()
plt.savefig("hierarchical_silhouette.png")
plt.close()

# Plot dendrogram
plt.figure(figsize=(12, 6))
dendrogram(Z, no_labels=True)
plt.title("Hierarchical Clustering Dendrogram")
plt.xlabel("Sample")
plt.ylabel("Distance")
plt.tight_layout()
plt.savefig("hierarchical_dendrogram.png")
plt.show()

plt.close()
print("Saved dendrogram to hierarchical_dendrogramtest.png")

# Plot clusters on PCA projection
plt.figure(figsize=(8, 6))
plt.scatter(X[:, 0], X[:, 1], c=best_labels, cmap='tab10', s=35, alpha=0.8)
plt.title(f"Hierarchical Clusters on PCA (k={best_k})")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.grid(True)
plt.colorbar(label="Cluster")
plt.tight_layout()
plt.savefig("hierarchical_pca_clusters.png")
plt.show()

plt.close()
print("Saved PCA cluster plot to hierarchical_pca_clusters.png")