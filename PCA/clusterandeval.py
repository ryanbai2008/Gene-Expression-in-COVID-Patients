import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.mixture import GaussianMixture
from scipy.cluster.hierarchy import linkage, fcluster
import seaborn as sns
import matplotlib.pyplot as plt

# ==============================
# LOAD DATA
# ==============================
df = pd.read_csv("./data/normalized_counts_DESeq2(2).csv", index_col=0)
metadata = pd.read_excel("data/data.xlsx", index_col=0)

# Transpose so rows = samples, columns = genes
X = df.T
X = X.apply(pd.to_numeric, errors='coerce').fillna(0)

# ==============================
# FEATURE SELECTION (Variance)
# ==============================
selector = VarianceThreshold(threshold=0.37)
X_highvar = selector.fit_transform(X)
selected_genes = X.columns[selector.get_support()]
X_highvar_df = pd.DataFrame(X_highvar, index=X.index, columns=selected_genes)

print(f"Number of genes after variance filtering: {X_highvar_df.shape[1]}")

# ==============================
# STANDARDIZE
# ==============================
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_highvar_df)

# ==============================
# PCA
# ==============================
pca = PCA(n_components=20)
X_pca = pca.fit_transform(X_scaled)

total_var = np.sum(pca.explained_variance_ratio_) * 100
print(f"Total variance captured by 20 PCs: {total_var:.2f}%")

# Scree plot to visualize cumulative variance
plt.figure(figsize=(6,4))
plt.plot(np.cumsum(pca.explained_variance_ratio_) * 100, marker='o')
plt.xlabel("Number of PCs")
plt.ylabel("Cumulative Explained Variance (%)")
plt.title("PCA Scree Plot")
plt.grid(True)
plt.tight_layout()
plt.show()

# ==============================
# CLUSTERING COMPARISON
# ==============================
methods = ["Ward", "Average", "KMeans", "GaussianMixture"]
num_clusters = 3
n_components = 5  # using first 5 PCs for clustering
X_pca_sub = X_pca[:, :n_components]

fig, axes = plt.subplots(1, len(methods), figsize=(20, 5), sharex=True, sharey=True)

for i, method in enumerate(methods):
    if method in ["Ward", "Average"]:
        linked = linkage(X_pca_sub, method=method.lower())
        cluster_labels = fcluster(linked, t=num_clusters, criterion='maxclust')
    elif method == "KMeans":
        cluster_labels = KMeans(n_clusters=num_clusters, random_state=42).fit_predict(X_pca_sub)
    elif method == "GaussianMixture":
        cluster_labels = GaussianMixture(n_components=num_clusters, random_state=42).fit_predict(X_pca_sub)

    pca_df = pd.DataFrame(X_pca_sub[:, :2], columns=["PC1", "PC2"], index=X.index)
    pca_df["Cluster"] = cluster_labels

    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Cluster", palette="Set2", ax=axes[i], s=60)
    axes[i].set_title(f"{method}")
    axes[i].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)")
    axes[i].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)")
    axes[i].legend().remove()

plt.suptitle(f"Comparison of Clustering Methods on PCA ({n_components} PCs)", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()

# ==============================
# CLUSTER QUALITY METRICS
# ==============================
cluster_methods = {
    "Ward (Hierarchical)": AgglomerativeClustering(n_clusters=num_clusters, linkage="ward"),
    "Average (Hierarchical)": AgglomerativeClustering(n_clusters=num_clusters, linkage="average"),
    "KMeans": KMeans(n_clusters=num_clusters, random_state=42),
    "GaussianMixture": GaussianMixture(n_components=num_clusters, random_state=42),
}

scores = []
for name, model in cluster_methods.items():
    labels = model.fit_predict(X_pca_sub)
    sil = silhouette_score(X_pca_sub, labels)
    ch = calinski_harabasz_score(X_pca_sub, labels)
    db = davies_bouldin_score(X_pca_sub, labels)
    scores.append([name, sil, ch, db])

score_df = pd.DataFrame(scores, columns=["Method", "Silhouette", "Calinski-Harabasz", "Davies-Bouldin"])
print("\n=== Cluster Evaluation Metrics ===")
print(score_df.sort_values("Silhouette", ascending=False).to_string(index=False))

# ==============================
# OPTIMAL NUMBER OF PCs TEST
# ==============================
scores = []
for n in range(2, 21):  # test 2–20 PCs
    X_pca_sub = X_pca[:, :n]
    labels = KMeans(n_clusters=3, random_state=42).fit_predict(X_pca_sub)
    sil = silhouette_score(X_pca_sub, labels)
    scores.append(sil)

plt.figure(figsize=(6,4))
plt.plot(range(2, 21), scores, marker='o')
plt.xlabel("Number of PCs used")
plt.ylabel("Silhouette Score")
plt.title("Optimal Number of PCs for Clustering (KMeans)")
plt.grid(True)
plt.tight_layout()
plt.show()
