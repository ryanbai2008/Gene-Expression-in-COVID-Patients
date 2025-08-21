import pandas as pd
from sklearn.cluster import SpectralClustering

data_df = pd.read_csv("data2/normalized_counts_DESeq2(2).csv", index_col=0)
data_df = data_df.T

spectral = SpectralClustering(n_clusters=5, affinity="nearest_neighbors", random_state=10, n_init=3)
kernel_kmeans_labels = spectral.fit_predict(data_df.values)

clusters_df = pd.DataFrame({"Patient" : data_df.index, "Cluster" : kernel_kmeans_labels})
clusters_df.to_csv("models/kernel_KMeans_clusters.csv", index=False)
