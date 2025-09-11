import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import seaborn as sns
import re

# -----------------------------
# Load expression data
# -----------------------------
df = pd.read_csv("./data/normalized_counts_DESeq2(2).csv", index_col=0)
X = df.T.copy()  # samples as rows

# Normalize sample names
def normalize_sample_name(name):
    return re.sub(r'd\d+$', '', name)
X['Normalized'] = X.index.to_series().apply(normalize_sample_name)

# Load metadata
metadata = pd.read_excel("data/data.xlsx")
metadata = metadata.set_index('Patient_ID')

# Merge metadata
X = X.merge(metadata, left_on='Normalized', right_index=True, how='inner')
group_labels = X['Status'].copy()

# Drop extra columns
X.drop(columns=['Normalized', 'Status'], inplace=True)
X = X.apply(pd.to_numeric, errors='coerce').fillna(0)

# -----------------------------
# Filter high-variance genes
# -----------------------------
selector = VarianceThreshold(threshold=0.37)
X_highvar = selector.fit_transform(X)
selected_genes = X.columns[selector.get_support()]
X_highvar_df = pd.DataFrame(X_highvar, index=X.index, columns=selected_genes)

# Standardize
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_highvar_df)

# -----------------------------
# Subset to unhealthy samples only
# -----------------------------
unhealthy_types = [
    'COVID-19',
    'Influenza_A',
    'Influenza_B',
    'Coronavirus',
    'Co-infection (Bac/viral)',
    'Co-infection (Viral/Fungal)',
    'Co-infection (Bac/Viral/Fungal)',
    'Co-infection (viral/viral/fungal)',
    'Co-infection (Viral/viral)',
    'Sepsis',
    'Septic_shock'
]

unhealthy_mask = group_labels.isin(unhealthy_types)
X_unhealthy = X_scaled[unhealthy_mask.values, :]
unhealthy_indices = X_highvar_df.index[unhealthy_mask]
unhealthy_labels = group_labels[unhealthy_mask]

print(f"Unhealthy subset shape: {X_unhealthy.shape}")

# -----------------------------
# Test multiple k and compute silhouette
# -----------------------------
K_range = range(2, 7)  # test 2 to 6 clusters
sil_scores = []


pca = PCA(n_components=10)
X_pca = pca.fit_transform(X_scaled)
print(f"Explained variance by PC1 and PC2: {pca.explained_variance_ratio_[:10].sum()*100:.2f}%")

for k in K_range:
    kmeans = KMeans(n_clusters=k, random_state=10, n_init=10)
    labels = kmeans.fit_predict(X_pca)
    score = silhouette_score(X_pca, labels)
    sil_scores.append(score)
    print(f"k={k}, silhouette score={score:.3f}")

# Plot silhouette vs k
plt.figure(figsize=(6,4))
plt.plot(K_range, sil_scores, marker='o')
plt.xlabel("Number of clusters (k)")
plt.ylabel("Average Silhouette Score")
plt.title("Silhouette Analysis for Unhealthy Samples")
plt.show()

# -----------------------------
# Choose optimal k
# -----------------------------
best_k = K_range[np.argmax(sil_scores)]
print(f"Best k based on silhouette: {best_k}")

# -----------------------------
# Run final K-means with best k
# -----------------------------
kmeans_final = KMeans(n_clusters=best_k, random_state=10, n_init=10)
unhealthy_labels_cluster = kmeans_final.fit_predict(X_unhealthy)

# Map labels back to full dataset
full_labels = pd.Series(data=np.nan, index=X_highvar_df.index)
full_labels[unhealthy_indices] = unhealthy_labels_cluster

# -----------------------------
# PCA for plotting
# -----------------------------
pca = PCA(n_components=50)
X_pca = pca.fit_transform(X_scaled)
pc1_var = pca.explained_variance_ratio_[0]*100
pc2_var = pca.explained_variance_ratio_[1]*100

# Plot PCA with all samples colored by infection type
unique_groups = group_labels.unique()
palette = sns.color_palette("tab10", n_colors=len(unique_groups))
group_color_map = dict(zip(unique_groups, palette))

plt.figure(figsize=(10,7))
for grp in unique_groups:
    idx = group_labels[group_labels == grp].index
    plt.scatter(X_pca[X_highvar_df.index.get_indexer(idx),0],
                X_pca[X_highvar_df.index.get_indexer(idx),1],
                label=grp,
                alpha=0.6,
                color=group_color_map[grp])

# Overlay clusters for unhealthy samples
for cl in np.unique(unhealthy_labels_cluster):
    cluster_idx = unhealthy_indices[unhealthy_labels_cluster == cl]
    plt.scatter(X_pca[X_highvar_df.index.get_indexer(cluster_idx),0],
                X_pca[X_highvar_df.index.get_indexer(cluster_idx),1],
                edgecolor='black', s=100, facecolor='none',
                linewidth=1.5,
                label=f"Unhealthy cluster {cl}")

plt.xlabel(f"PC1 ({pc1_var:.1f}% variance)")
plt.ylabel(f"PC2 ({pc2_var:.1f}% variance)")
plt.title("PCA of Samples with Unhealthy Clusters Highlighted")
plt.legend(bbox_to_anchor=(1.05,1), loc='upper left')
plt.tight_layout()
plt.show()

# -----------------------------
# Save cluster assignments
# -----------------------------
cluster_df = pd.DataFrame({
    "Sample": unhealthy_indices,
    "Cluster": unhealthy_labels_cluster,
    "Infection_Status": unhealthy_labels.values
})
cluster_df.to_csv("unhealthy_clusters_bestk.csv", index=False)
