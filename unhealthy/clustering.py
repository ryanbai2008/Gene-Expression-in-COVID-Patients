import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold
import matplotlib.pyplot as plt
import seaborn as sns
import re
from numpy.linalg import norm

# -----------------------------
# Load and prepare expression data
# -----------------------------
df = pd.read_csv("./data/normalized_counts_DESeq2(2).csv", index_col=0)
X = df.T.copy()  # rows = samples, cols = genes

# Normalize sample names to match metadata
def normalize_sample_name(name):
    return re.sub(r'd\d+$', '', name)

X['Normalized'] = X.index.to_series().apply(normalize_sample_name)

# Load metadata
metadata = pd.read_excel("data/data.xlsx")
metadata = metadata.set_index('Patient_ID')

# Merge metadata with expression data
X = X.merge(metadata, left_on='Normalized', right_index=True, how='inner')
group_labels = X['Status'].copy()  # infection status: COVID, flu, co-infection, etc.

# Drop extra columns
X.drop(columns=['Normalized', 'Status'], inplace=True)

# Ensure numeric and fill NaNs
X = X.apply(pd.to_numeric, errors='coerce').fillna(0)

# -----------------------------
# Filter high-variance genes
# -----------------------------
selector = VarianceThreshold(threshold=0.37)
X_highvar = selector.fit_transform(X)
selected_genes = X.columns[selector.get_support()]
X_highvar_df = pd.DataFrame(X_highvar, index=X.index, columns=selected_genes)

# -----------------------------
# Standardize features
# -----------------------------
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_highvar_df)

# -----------------------------
# PCA on all samples (for reference plots)
# -----------------------------
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)
pc1_var = pca.explained_variance_ratio_[0] * 100
pc2_var = pca.explained_variance_ratio_[1] * 100

# -----------------------------
# Role 12: Unhealthy-only clustering
# -----------------------------
unhealthy_statuses = [
    'COVID-19', 'Influenza_A', 'Influenza_B', 
    'Co-infection (Bac/viral)', 'Co-infection (Viral/Fungal)',
    'Co-infection (Bac/Viral/Fungal)', 'Co-infection (viral/viral/fungal)', 'Co-infection (Viral/viral)', 'Sepsis', 'Septic_shock',
    'Coronavirus'
]

unhealthy_mask = group_labels.isin(unhealthy_statuses)
X_unhealthy = X_scaled[unhealthy_mask.values, :]
unhealthy_indices = X_highvar_df.index[unhealthy_mask]

# Run K-means
kmeans_unhealthy = KMeans(n_clusters=2, random_state=10, n_init=10)
unhealthy_labels = kmeans_unhealthy.fit_predict(X_unhealthy)

# -----------------------------
# Role 13: Outlier detection
# -----------------------------
centroids = np.array([X_unhealthy[unhealthy_labels==i].mean(axis=0) 
                      for i in np.unique(unhealthy_labels)])

distances = np.array([norm(X_unhealthy[i]-centroids[unhealthy_labels[i]]) 
                      for i in range(X_unhealthy.shape[0])])

threshold = distances.mean() + 2*distances.std()
outliers = unhealthy_indices[distances > threshold]
print("Unhealthy outlier samples:", outliers.tolist())

# -----------------------------
# Role 14: Map clusters to infection status
# -----------------------------
cluster_df = pd.DataFrame({
    "Sample": unhealthy_indices,
    "Cluster": unhealthy_labels,
    "Infection": group_labels[unhealthy_mask].values
})

# Optional: contingency table
contingency = pd.crosstab(cluster_df['Cluster'], cluster_df['Infection'])
print("Cluster vs Infection status:\n", contingency)

# -----------------------------
# PCA plot: unhealthy clusters
# -----------------------------
plt.figure(figsize=(8,6))
unique_clusters = np.unique(unhealthy_labels)
palette = sns.color_palette("tab10", n_colors=len(unique_clusters))

for cl, color in zip(unique_clusters, palette):
    idx = unhealthy_indices[unhealthy_labels == cl]
    plt.scatter(X_pca[X_highvar_df.index.get_indexer(idx), 0],
                X_pca[X_highvar_df.index.get_indexer(idx), 1],
                label=f"Cluster {cl}", alpha=0.8, color=color)

plt.xlabel(f"PC1 ({pc1_var:.1f}% variance)")
plt.ylabel(f"PC2 ({pc2_var:.1f}% variance)")
plt.title("PCA of Unhealthy Samples (Clusters)")
plt.legend()
plt.tight_layout()
plt.show()

# -----------------------------
# Optional: Save cluster assignments
# -----------------------------
cluster_df.to_csv("unhealthy_clusters.csv", index=False)
