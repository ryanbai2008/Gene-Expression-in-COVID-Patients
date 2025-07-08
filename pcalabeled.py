import pandas as pd
import numpy as np
import re
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold
import matplotlib.pyplot as plt
import seaborn as sns

# Load expression matrix (genes as index, samples as columns)
# DO EITHER 2 or 4
df = pd.read_csv("./data2/normalized_counts_DESeq2(2).csv", index_col=0)

# Transpose so rows = samples, columns = genes
X = df.T.copy()

# Normalize sample names (e.g., nCWES001d1 → nCWES001)
def normalize_sample_name(name):
    return re.sub(r'd\d+$', '', name)

print(X[150:220])
X['Normalized'] = X.index.to_series().apply(normalize_sample_name)
print(X[150:220])

# Load metadata
metadata = pd.read_excel("data/data.xlsx")
metadata = metadata.set_index('Patient_ID')  # Ensure 'Patient_ID' is the index
print(metadata['Status'].unique())

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

# Plot PCA scatter with labels
unique_groups = group_labels.unique()
palette = sns.color_palette("tab10", n_colors=len(unique_groups))
group_color_map = dict(zip(unique_groups, palette))

plt.figure(figsize=(8, 6))
for grp in unique_groups:
    idx = group_labels[group_labels == grp].index
    plt.scatter(X_pca[X_highvar_df.index.get_indexer(idx), 0],
                X_pca[X_highvar_df.index.get_indexer(idx), 1],
                label=grp, alpha=0.8, color=group_color_map[grp])

pc1_var = pca.explained_variance_ratio_[0] * 100
pc2_var = pca.explained_variance_ratio_[1] * 100

plt.xlabel(f"PC1 ({pc1_var:.1f}% variance)")
plt.ylabel(f"PC2 ({pc2_var:.1f}% variance)")
plt.title("PCA of High-Variance Genes")
plt.legend(title="Condition")
plt.tight_layout()
plt.show()

# Print variance explained and number of genes
print(f"Explained variance by PC1 and PC2: {pca.explained_variance_ratio_[:2].sum()*100:.2f}%")
print(f"Number of genes after variance filtering: {X_highvar_df.shape[1]}")
