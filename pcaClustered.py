import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.feature_selection import VarianceThreshold
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import fcluster


# Load expression matrix (genes as index, samples as columns)
df = pd.read_csv("./data2/normalized_counts_DESeq2(2).csv", index_col=0)
metadata = pd.read_excel("data/data.xlsx", index_col=0)  

# Transpose so rows = samples, columns = genes (PCA expects features per sample)
X = df.T

print(X.index[:5])
print(metadata.index[:5])

print(X.dtypes.value_counts())
X = X.apply(pd.to_numeric, errors='coerce')  # Converts strings to NaN
X = X.fillna(0)  # Replace any remaining NaNs with 0

selector = VarianceThreshold(threshold=0.37)  
X_highvar = selector.fit_transform(X)

# Get the names of the selected genes
selected_genes = X.columns[selector.get_support()]
X_highvar_df = pd.DataFrame(X_highvar, index=X.index, columns=selected_genes)

# Standardize features (mean=0, std=1)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_highvar_df)

pca = PCA(n_components=2)  
X_pca = pca.fit_transform(X_scaled)

#to cluster
linked = linkage(X_pca, method='ward')  # or use 'average', 'complete', etc.

# Plot dendrogram
# plt.figure(figsize=(12, 6))
# dendrogram(linked, labels=X.index.to_list(), leaf_rotation=90)
# plt.title('Hierarchical Clustering Dendrogram (on PCA)')
# plt.xlabel('Sample')
# plt.ylabel('Distance')
# plt.tight_layout()
# plt.show()


# Print explained variance
print(f"Explained variance by PC1 and PC2: {pca.explained_variance_ratio_[:2].sum()*100:.2f}%")
print(f"Number of genes after variance filtering: {X_highvar_df.shape[1]}")

#clustering data
cluster_labels = fcluster(linked, t=2, criterion='maxclust')


# Create dataframe for plotting
pca_df = pd.DataFrame(X_pca, columns=['PC1', 'PC2'], index=X.index)
pca_df['Sample'] = X.index
pca_df['Cluster'] = cluster_labels



# Plot results
plt.figure(figsize=(8,6))
sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Cluster', palette='Set2')
plt.title('PCA of GSE217948 Expression Matrix with Clustering')
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)")
plt.grid(True)
plt.tight_layout()
plt.show()
