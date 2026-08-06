# prevent memory leak issues n stuff
import os
os.environ['OMP_NUM_THREADS'] = '1'

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# get data, make sure to transpose bc we need it in that format
data_df = pd.read_csv("normalized_gene_data_log2.csv", index_col=0)
data_df = data_df.T
print(data_df.shape)
# data_df = data_df.iloc[1:, :]

# scale the data before clustering
scaler = StandardScaler()
X_scaled = scaler.fit_transform(data_df)


pca = PCA(n_components=10, random_state=42)
X_pca = pca.fit_transform(X_scaled)

explained_variance = pca.explained_variance_ratio_
cumulative_variance = np.cumsum(explained_variance)

# fit KMeans
k = 2
kmeans = KMeans(n_clusters=k, random_state=42, n_init=3)
kmeans_labels = kmeans.fit_predict(X_pca)

print(kmeans_labels.shape)
print(kmeans_labels)

# save cluster assignments
cluster_df = pd.DataFrame({'cluster': kmeans_labels}, index=data_df.index)
cluster_file = os.path.join(os.path.dirname(__file__), f'kmeans_clusters_k{k}.csv')
cluster_df.to_csv(cluster_file)
print(f'Saved KMeans cluster assignments to {cluster_file}')

# PCA for 2D visualization
pca = PCA(n_components=2, random_state=42)
X_pca2 = pca.fit_transform(X_pca)

plt.figure(figsize=(8, 6))
plt.scatter(X_pca2[:, 0], X_pca2[:, 1], c=kmeans_labels, cmap='viridis', s=35)
plt.title('K-Means Clusters (PCA projection)')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.colorbar(label='Cluster')
plt.tight_layout()
plt.savefig(os.path.join(os.path.dirname(__file__), 'kmeans_pca_plotk=2.png'))
plt.show()

# elbow plot for a quick k selection view
wcss = []
for cluster_count in range(1, 11):
    model = KMeans(n_clusters=cluster_count, random_state=10, n_init=3)
    model.fit(X_scaled)
    wcss.append(model.inertia_)

plt.figure(figsize=(8, 5))
plt.plot(range(1, 11), wcss, marker='o')
plt.title('Elbow Method for K-Means')
plt.xlabel('Number of clusters (k)')
plt.ylabel('Within-Cluster Sum of Squares (WCSS)')
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(os.path.dirname(__file__), 'kmeans_elbow_plot1.png'))
plt.show()
