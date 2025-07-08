import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.feature_selection import VarianceThreshold

# Load expression matrix (genes as index, samples as columns)
df = pd.read_csv("data/GSE217948_combined_expression_matrix.csv")

# Transpose so rows = samples, columns = genes (PCA expects features per sample)
X = df.T

print(X.dtypes.value_counts())
X = X.apply(pd.to_numeric, errors='coerce')  # Converts strings to NaN
X = X.fillna(0)  # Replace any remaining NaNs with 0

selector = VarianceThreshold(threshold=20)  
X_highvar = selector.fit_transform(X)

# Get the names of the selected genes
selected_genes = X.columns[selector.get_support()]
X_highvar_df = pd.DataFrame(X_highvar, index=X.index, columns=selected_genes)

# Log-transform (log1p avoids issues with 0s)
X_log = np.log1p(X_highvar_df)


# Standardize features (mean=0, std=1)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_log)

pca = PCA(n_components=2)  
X_pca = pca.fit_transform(X_scaled)

# Print explained variance
print(f"Explained variance by PC1 and PC2: {pca.explained_variance_ratio_[:2].sum()*100:.2f}%")

# Create dataframe for plotting
pca_df = pd.DataFrame(X_pca, columns=['PC1', 'PC2'])
pca_df['Sample'] = X.index

# Plot results
plt.figure(figsize=(8,6))
sns.scatterplot(data=pca_df, x='PC1', y='PC2')
plt.title('PCA of GSE217948 Expression Matrix')
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)")
plt.grid(True)
plt.tight_layout()
plt.show()
