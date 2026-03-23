from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
import seaborn as sns
import re

# Load expression matrix
df = pd.read_csv("./data2/normalized_counts_DESeq2(2).csv", index_col=0)
X = df.T.copy()

# Normalize sample names
def normalize_sample_name(name):
    return re.sub(r'd\d+$', '', name)

X['Normalized'] = X.index.to_series().apply(normalize_sample_name)

# Load metadata and merge
metadata = pd.read_excel("data/data.xlsx").set_index('Patient_ID')
X = X.merge(metadata, left_on='Normalized', right_index=True, how='inner')
group_labels = X['Status'].copy()
X.drop(columns=['Normalized', 'Status'], inplace=True)

# Convert to numeric and fill NaNs
X = X.apply(pd.to_numeric, errors='coerce').fillna(0)

# High variance gene filtering
selector = VarianceThreshold(threshold=0.37)
X_highvar = selector.fit_transform(X)
selected_genes = X.columns[selector.get_support()]
X_highvar_df = pd.DataFrame(X_highvar, index=X.index, columns=selected_genes)

# Standardize
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_highvar_df)

# Define function to run t-SNE and score
def run_tsne(data, use_pca=False, n_pca_components=50, title=''):
    if use_pca:
        pca = PCA(n_components=n_pca_components)
        data = pca.fit_transform(data)
    tsne = TSNE(n_components=2, perplexity=20, max_iter=1000, random_state=42)
    embedding = tsne.fit_transform(data)
    
    # KMeans to get labels for silhouette scoring
    kmeans = KMeans(n_clusters=len(group_labels.unique()), random_state=42)
    km_labels = kmeans.fit_predict(embedding)
    score = silhouette_score(embedding, km_labels)
    
    return embedding, score

# Run both
tsne_raw, score_raw = run_tsne(X_scaled, use_pca=False, title="t-SNE (No PCA)")
tsne_pca, score_pca = run_tsne(X_scaled, use_pca=True, n_pca_components=50, title="t-SNE (PCA=50)")

# Plot both
palette = sns.color_palette("tab10", n_colors=len(group_labels.unique()))
group_color_map = dict(zip(group_labels.unique(), palette))

fig, axs = plt.subplots(1, 2, figsize=(16, 7))

for i, (tsne, score, title, ax) in enumerate(zip(
    [tsne_raw, tsne_pca],
    [score_raw, score_pca],
    ["t-SNE (No PCA)", "t-SNE (PCA=50)"],
    axs
)):
    for grp in group_labels.unique():
        idx = group_labels[group_labels == grp].index
        locs = X_highvar_df.index.get_indexer(idx)
        ax.scatter(tsne[locs, 0], tsne[locs, 1], label=grp, alpha=0.8, color=group_color_map[grp])
    ax.set_title(f"{title}\nSilhouette Score = {score:.3f}")
    ax.set_xlabel("t-SNE 1")
    ax.set_ylabel("t-SNE 2")
    ax.legend(loc='best', fontsize='small')

plt.tight_layout()
plt.show()
