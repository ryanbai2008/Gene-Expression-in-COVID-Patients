import pandas as pd
import numpy as np
from sklearn.feature_selection import VarianceThreshold

data = pd.read_csv('data2/normalized_counts_DESeq2(2).csv', index_col=0)
data = data.iloc[:-4]
data = data.T

selector = VarianceThreshold(threshold=0.37)
X_highvar = selector.fit_transform(data)
selected_genes = data.columns[selector.get_support()]
X_df = pd.DataFrame(X_highvar, index=data.index, columns=selected_genes)

def pearson_correlation(data):
    standardized = (data - np.mean(data, axis=0)) / np.std(data, axis=0)
    return np.corrcoef(standardized, rowvar=False)

def adjacency(cor_matrix, beta):
    return np.abs(cor_matrix) ** beta

def compute_tom(adjacency_matrix):
    k = np.sum(adjacency_matrix, axis=1)
    num_genes = adjacency_matrix.shape[0]
    shared = adjacency_matrix @ adjacency_matrix
    min_k = np.minimum.outer(k, k)
    np.fill_diagonal(min_k, k)
    denom = min_k - adjacency_matrix + 1e-10
    tom = (shared + adjacency_matrix) / denom
    np.fill_diagonal(tom, 1.0)
    return tom

cor = pearson_correlation(X_df.to_numpy())
adj = adjacency(cor, beta=6)
tom = 1 - compute_tom(adj)
tom = np.clip(tom, 0, 1)

import json
tom_df = pd.DataFrame(tom, index=selected_genes, columns=selected_genes)
corr_matrix = pd.DataFrame(cor, index=selected_genes, columns=selected_genes)
corr_matrix.to_csv("correlation_matrix_values.csv")
tom_df.to_csv("tom_dissimilarity_values.csv")
with open("Genes used for analysis (TOM).json", 'w') as file:
    json.dump(selected_genes.tolist(), file)
