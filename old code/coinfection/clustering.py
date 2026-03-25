import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import re

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
from sklearn.feature_selection import VarianceThreshold
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
# PCA on all samples (full dataset)
# -----------------------------
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)
pc1_var = pca.explained_variance_ratio_[0] * 100
pc2_var = pca.explained_variance_ratio_[1] * 100

# -----------------------------
# Run k-means on co-infection subset (Role 15 output)
# -----------------------------
co_mask = group_labels.str.contains("Co-infection", case=False, na=False)
X_co = X_scaled[co_mask.values, :]
co_indices = X_highvar_df.index[co_mask]
kmeans_co = KMeans(n_clusters=3, random_state=10, n_init=10)
co_labels = kmeans_co.fit_predict(X_co)

# Map co-infection cluster labels back to full dataset
full_labels = pd.Series(data=np.nan, index=X_highvar_df.index)
full_labels[co_indices] = co_labels

# -----------------------------
# Plot PCA: all samples colored by infection type
# Overlay co-infection clusters as special markers
# -----------------------------
unique_groups = group_labels.unique()
palette = sns.color_palette("tab10", n_colors=len(unique_groups))
group_color_map = dict(zip(unique_groups, palette))

plt.figure(figsize=(8,6))
for grp in unique_groups:
    idx = group_labels[group_labels == grp].index
    plt.scatter(X_pca[X_highvar_df.index.get_indexer(idx), 0],
                X_pca[X_highvar_df.index.get_indexer(idx), 1],
                label=grp,
                alpha=0.6,
                color=group_color_map[grp])

# Overlay co-infection clusters
for cl in np.unique(co_labels):
    cluster_idx = co_indices[co_labels == cl]
    plt.scatter(X_pca[X_highvar_df.index.get_indexer(cluster_idx), 0],
                X_pca[X_highvar_df.index.get_indexer(cluster_idx), 1],
                edgecolor='black', s=100, facecolor='none',
                linewidth=1.5,
                label=f"Co-infection cluster {cl}")

plt.xlabel(f"PC1 ({pc1_var:.1f}% variance)")
plt.ylabel(f"PC2 ({pc2_var:.1f}% variance)")
plt.title("PCA of Samples with Co-infection Clusters Highlighted")
plt.legend(bbox_to_anchor=(1.05,1), loc='upper left')
plt.tight_layout()
plt.show()

# -----------------------------
# Optional: save co-infection cluster assignments
# -----------------------------
co_cluster_df = pd.DataFrame({
    "Sample": co_indices,
    "Cluster": co_labels
})
co_cluster_df.to_csv("co_infection_clusters_fullPCA.csv", index=False)
