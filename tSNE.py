from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.datasets import make_classification
import numpy as np
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import seaborn as sns
from sklearn.feature_selection import VarianceThreshold
import re


# Generate some sample data

df = pd.read_csv("./data2/normalized_counts_DESeq2(2).csv", index_col=0)

# Transpose so rows = samples, columns = genes
X = df.T.copy()

# Normalize sample names (e.g., nCWES001d1 → nCWES001)
def normalize_sample_name(name):
    return re.sub(r'd\d+$', '', name)

X['Normalized'] = X.index.to_series().apply(normalize_sample_name)

# Load metadata
metadata = pd.read_excel("data/data.xlsx")
metadata = metadata.set_index('Patient_ID')  

# Merge metadata based on normalized names
X = X.merge(metadata, left_on='Normalized', right_index=True, how='inner')
group_labels = X['Status'].copy()

# Drop extra columns (Normalized + Status)
X.drop(columns=['Normalized', 'Status'], inplace=True)

# Ensure numeric, fill NaNs
X = X.apply(pd.to_numeric, errors='coerce').fillna(0)

# Filter for high variance genes
selector = VarianceThreshold(threshold=0.37)
X_highvar = selector.fit_transform(X)
selected_genes = X.columns[selector.get_support()]
X_highvar_df = pd.DataFrame(X_highvar, index=X.index, columns=selected_genes)

# Standardize features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_highvar_df)

# Perform PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)
# Initialize t-SNE with desired parameters
X_tsne = TSNE(n_components=2, perplexity=20, max_iter=1000).fit_transform(X_pca) #SAM CHANGE THESE PARAMETERS for (Perplexity, max_iter) for better clustering

# Fit t-SNE to the data
unique_groups = group_labels.unique()
palette = sns.color_palette("tab10", n_colors=len(unique_groups))
group_color_map = dict(zip(unique_groups, palette))

plt.figure(figsize=(8,6))

for grp in unique_groups:
    idx = group_labels[group_labels == grp].index
    locs = X_highvar_df.index.get_indexer(idx)
    plt.scatter(X_tsne[locs, 0], X_tsne[locs, 1], label=grp, alpha=0.8, color=group_color_map[grp])

plt.title("t-SNE colored by Status")
plt.xlabel("t-SNE 1")
plt.ylabel("t-SNE 2")
plt.legend(title="Status")
plt.show()
