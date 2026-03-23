import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
df = pd.read_csv("data/normalized_logCPM.csv")

# Transpose so PCA sees rows = samples
df = df.T

# Standardize the data (important for PCA)
scaler = StandardScaler()
scaled = scaler.fit_transform(df)

# Run PCA
pca = PCA(n_components=2)
components = pca.fit_transform(scaled)

# Create DataFrame for plotting
pca_df = pd.DataFrame(components, columns=['PC1', 'PC2'], index=df.index)

# Plot
plt.figure(figsize=(8,6))
sns.scatterplot(data=pca_df, x='PC1', y='PC2')
plt.title('PCA of Gene Expression')
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
plt.grid(True)
plt.tight_layout()
plt.show()
