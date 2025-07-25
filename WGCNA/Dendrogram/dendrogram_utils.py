from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
import pandas as pd
import numpy as np

TOM_dissimilarity_df = pd.read_csv("tom_dissimilarity_values.csv", index_col=0)
gene_names_list = TOM_dissimilarity_df.columns.tolist()
TOM_dissimilarity = TOM_dissimilarity_df.values
TOM_dissimilarity = (TOM_dissimilarity + TOM_dissimilarity.T) / 2
np.fill_diagonal(TOM_dissimilarity, 0)

upper_TOM_dissimilarity = squareform(TOM_dissimilarity, checks=True)

linkage_matrix = linkage(upper_TOM_dissimilarity, method='average')

threshold = 750
clusters = fcluster(linkage_matrix, t=threshold, criterion='maxclust')

# save the data to files
gene_clusters_df = pd.DataFrame({
    'Gene': TOM_dissimilarity_df.index,
    'Cluster': clusters
})

gene_clusters_df.sort_values(by='Cluster', inplace=True)
gene_clusters_df.to_csv("gene_clusters.csv", index=False)

