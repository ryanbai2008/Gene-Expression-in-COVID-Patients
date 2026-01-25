import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import sys
import os

os.environ["OMP_NUM_THREADS"] = "2"

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from data.getData import *

# data!!
X, y = getData()

# scaling
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# setting k (from elbow method used previously)
k = 3

# fitting kmeans
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)
kmeans = KMeans(n_clusters=k, random_state=42, n_init='auto')
clusters = kmeans.fit_predict(X_scaled)

print("cluster labels for each sample:")
print(clusters)


plt.figure(figsize=(8,6))
scatter = plt.scatter(X_pca[:,0], X_pca[:,1], c=clusters, cmap='tab10', s=10)
plt.title('KMeans clustering w/k = 3')
plt.xlabel('PCA component 1')
plt.ylabel('PCA component 2')
plt.colorbar(scatter, label='Cluster ID')
plt.grid(True)
plt.tight_layout()
plt.show()
