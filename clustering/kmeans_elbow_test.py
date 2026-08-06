import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

os.environ["OMP_NUM_THREADS"] = "2"

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA

# data!!
data_df = pd.read_csv("normalized_gene_data_log2.csv", index_col=0)
data_df = data_df.T

print(f"Data shape: {data_df.shape}")

# HVG Selection (not using since it is better to use full genes)
'''
gene_mean = data_df.mean(axis=0)
gene_var = data_df.var(axis=0)

bins = pd.qcut(gene_mean, q=20, duplicates='drop')
stats = pd.DataFrame({'mean': gene_mean, 'var': gene_var, 'bin': bins})
stats['var_z'] = stats.groupby('bin', observed=True)['var'].transform(
    lambda x: (x - x.mean()) / (x.std(ddof=0) + 1e-12)
)

top_genes = stats['var_z'].nlargest(2000).index

data_df = data_df[top_genes]

print(f"Data shape: {data_df.shape}")
'''

# scaling
scaler = StandardScaler()
X_scaled = scaler.fit_transform(data_df)

pca = PCA(n_components=10, random_state=42)
X_pca = pca.fit_transform(X_scaled)

explained_variance = pca.explained_variance_ratio_
cumulative_variance = np.cumsum(explained_variance)

print("PCA explained variance per component:")
for pc_idx, variance in enumerate(explained_variance, start=1):
    print(f"  PC{pc_idx}: {variance:.4%}")

print("Cumulative explained variance:")
for pc_idx, variance in enumerate(cumulative_variance, start=1):
    print(f"  PC{pc_idx}: {variance:.4%}")

print(f"Data shape: {data_df.shape}")


# range of k values
k_vals = range(1, 11)
wcss = []
silhouette_scores = []

# find wcss and silhouette scores
for k in k_vals:
    kmeans = KMeans(n_clusters=k, random_state=42, n_init='auto')
    kmeans.fit(X_pca)
    wcss.append(kmeans.inertia_)

    if k >= 2:
        silhouette_scores.append(silhouette_score(X_pca, kmeans.labels_, metric='euclidean'))

# getting elbow
def findElbow(wcss):
    k_vals = np.arange(1, len(wcss) + 1)
    lineStart = np.array([k_vals[0], wcss[0]])
    lineEnd = np.array([k_vals[-1], wcss[-1]])

    dist = []
    for i in range(len(k_vals)):
        pt = np.array([k_vals[i], wcss[i]])
        # dist from pt to line
        num = np.abs(np.cross(lineEnd - lineStart, lineStart - pt))
        den = np.linalg.norm(lineEnd - lineStart)
        dist.append(num / den)

    elbowK = k_vals[np.argmax(dist)]
    return elbowK

elbow_k = findElbow(wcss)
silhouette_k_vals = np.arange(2, len(silhouette_scores) + 2)
best_silhouette_k = silhouette_k_vals[np.argmax(silhouette_scores)]

print(f"Elbow k: {elbow_k}")
print(f"Best k by silhouette score: {best_silhouette_k}")
print(f"Silhouette scores: {silhouette_scores}")

# plot wcss
plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(k_vals, wcss, marker='o')
plt.scatter(elbow_k, wcss[elbow_k - 1], color='green', s=80, label=f'Elbow at k={elbow_k}')
plt.title('Elbow Method for K-Means')
plt.xlabel('Number of clusters')
plt.ylabel('WCSS')
plt.grid(True)
plt.legend()

# plot silhouette scores
plt.subplot(2, 1, 2)
plt.plot(silhouette_k_vals, silhouette_scores, marker='o', color='purple')
plt.scatter(best_silhouette_k, silhouette_scores[np.argmax(silhouette_scores)], color='red', s=80,
            label=f'Best silhouette at k={best_silhouette_k}')
plt.title('Silhouette Scores for K-Means')
plt.xlabel('Number of clusters')
plt.ylabel('Silhouette score')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()

pc_numbers = np.arange(1, len(explained_variance) + 1)

# scree plot for variance
plt.figure(figsize=(8, 4))
plt.plot(pc_numbers, explained_variance, marker='o', linewidth=2, color='tab:blue')
plt.title('PCA Scree Plot')
plt.xlabel('Principal Component')
plt.ylabel('Explained Variance Ratio')
plt.grid(True)
plt.tight_layout()

plt.show()

