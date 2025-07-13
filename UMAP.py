import umap
import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets import load_digits
from sklearn.preprocessing import StandardScaler

# Load a sample dataset (e.g., the digits dataset)
digits = load_digits()
data = digits.data
target = digits.target

# Scale the data (important for UMAP)
scaler = StandardScaler()
data_scaled = scaler.fit_transform(data)

# Create a UMAP object and fit it to the data
reducer = umap.UMAP(random_state=42)  # You can adjust parameters like n_neighbors, min_dist, etc.
embedding = reducer.fit_transform(data_scaled)

# Plot the UMAP embedding
plt.figure(figsize=(8, 6))
plt.scatter(embedding[:, 0], embedding[:, 1], c=target, cmap='Spectral', s=5)
plt.gca().set_aspect('equal', 'datalim')
plt.colorbar(boundaries=np.arange(11)-0.5).set_ticks(np.arange(10))
plt.title('UMAP projection of the digits dataset', fontsize=12)
plt.show()
