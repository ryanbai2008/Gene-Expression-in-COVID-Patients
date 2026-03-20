import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.feature_selection import VarianceThreshold
from scipy.cluster.hierarchy import linkage, fcluster

# ==============================
# LOAD DATA
# ==============================
df = pd.read_csv("./data/normalized_counts_DESeq2(2).csv", index_col=0)
metadata = pd.read_excel("data/data.xlsx", index_col=0)

# Transpose so rows = samples, columns = genes
X = df.T

# Data cleaning
X = X.apply(pd.to_numeric, errors='coerce').fillna(0)

# ==============================
# FEATURE SELECTION (Variance)
# ==============================
selector = VarianceThreshold(threshold=0.37)
X_highvar = selector.fit_transform(X)
selected_genes = X.columns[selector.get_support()]
X_highvar_df = pd.DataFrame(X_highvar, index=X.index, columns=selected_genes)

# ==============================
# STANDARDIZE
# ==============================
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_highvar_df)

# ==============================
# PCA
# ==============================
pca = PCA(n_components=5)  # Use more components for clustering
X_pca = pca.fit_transform(X_scaled)

explained = np.sum(pca.explained_variance_ratio_[:2]) * 100
print(f"Explained variance by PC1–2: {explained:.2f}%")
print(f"Number of genes after variance filtering: {X_highvar_df.shape[1]}")

# ==============================
# CLUSTERING METHODS TO COMPARE
# ==============================
methods = ['ward', 'average', 'complete', 'single']
num_clusters = 3

fig, axes = plt.subplots(1, len(methods), figsize=(20, 5), sharex=True, sharey=True)

for i, method in enumerate(methods):
    linked = linkage(X_pca, method=method)
    cluster_labels = fcluster(linked, t=num_clusters, criterion='maxclust')

    pca_df = pd.DataFrame(X_pca, columns=[f'PC{i+1}' for i in range(X_pca.shape[1])], index=X.index)
    pca_df['Cluster'] = cluster_labels
    pca_df['Sample'] = X.index

    sns.scatterplot(
        data=pca_df, x='PC1', y='PC2', hue='Cluster',
        palette='Set2', ax=axes[i], s=60
    )
    axes[i].set_title(f"{method.capitalize()} linkage")
    axes[i].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)")
    axes[i].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)")
    axes[i].legend().remove()

plt.suptitle('Comparison of Hierarchical Clustering Methods on PCA (5 PCs)', fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()
