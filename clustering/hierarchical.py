import os
os.environ['OMP_NUM_THREADS'] = '1'

import pandas as pd
import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# Load your normalized gene expression data.
data_df = pd.read_csv("normalized_gene_data_log2.csv", index_col=False)
data_df = data_df.T

data_df = data_df.iloc[1:, :]

# Standardize before clustering.
scaler = StandardScaler()
X = scaler.fit_transform(data_df)



pca = PCA(n_components=10, random_state=42)
X = pca.fit_transform(X)

# Compute linkage matrix using Ward's method.
linkage_matrix = hierarchy.linkage(X, method="ward")

# Plot dendrogram.
plt.figure(figsize=(12, 6))
hierarchy.dendrogram(linkage_matrix, labels=data_df.index.tolist(), leaf_rotation=90)
plt.title("Hierarchical Clustering Dendrogram")
plt.xlabel("Sample")
plt.ylabel("Distance")
plt.tight_layout()
plt.savefig(os.path.join(os.path.dirname(__file__), "hierarchical_dendrogram.png"))
plt.show()

# Assign cluster labels from the linkage tree.
# Change t or criterion to pick the desired number of clusters.
cluster_labels = hierarchy.fcluster(linkage_matrix, t=2, criterion="maxclust")

result_df = pd.DataFrame({"sample": data_df.index, "cluster": cluster_labels})
result_df.to_csv(os.path.join(os.path.dirname(__file__), "hierarchical_clusters.csv"), index=False)
print(result_df.head())
print("Saved hierarchical cluster labels and dendrogram.")
