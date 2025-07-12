import numpy as np

def compute_tom(transformed_data, beta=6):
    pearson_correlation_matrix = pearson_correlation(transformed_data)
    adjacency_matrix = adjacency(pearson_correlation_matrix, beta)
    connectivity = np.sum(adjacency_matrix, axis=1)
    shared_neighbors = adjacency_matrix @ adjacency_matrix
    num_genes = adjacency_matrix.shape[0]
    TOM_matrix = np.zeros_like(adjacency_matrix)
    min_connectivity = np.minimum.outer(connectivity, connectivity)
    for i in range(num_genes):
        for j in range(num_genes):
            if i == j:
                TOM_matrix[i][j] = 1
            else:
                numerator = adjacency_matrix[i, j] + shared_neighbors[i, j]   
                denominator = 1 - adjacency_matrix[i, j] + min_connectivity[i, j]
                if denominator == 0:
                    TOM_matrix[i][j] = 0
                else:
                    TOM_matrix[i][j] = numerator / denominator
    return TOM_matrix

def pearson_correlation(transformed_data):
    pearson_correlation_matrix = transformed_data @ transformed_data.T
    return pearson_correlation_matrix

def adjacency(pearson_correlation_matrix, beta):
    return np.abs(pearson_correlation_matrix) ** beta



#calculating the actual TOM matrix, only need to do once as we will save the matrix in a file
import pandas as pd
transformed_data = pd.read_csv('data/normalized_counts.csv')
gene_names = transformed_data.iloc[:, 0]
transformed_data = transformed_data.iloc[:, 1:].to_numpy()
TOM_matrix = compute_tom(transformed_data)
tom_data_frame = pd.DataFrame(TOM_matrix, columns=gene_names, index=gene_names)
tom_data_frame.to_csv('tom_matrix_values.csv')
