import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.feature_selection import VarianceThreshold

# Load your dataset
dataset = pd.read_excel("data/data.xlsx")  # or pd.read_csv("yourfile.csv")

# Convert 'NA', 'no record', etc. to np.nan
dataset.replace(['NA', 'no record', 'Discharged alive'], np.nan, inplace=True)

# Drop columns that are text or metadata (not helpful for PCA directly)
drop_cols = ['Patient_ID', 'Gender', 'RecruitmentSite', 'location_recruit', 
             'location_recruit_o', 'Clinical_outcome', 'outcome']
dataset.drop(columns=drop_cols, inplace=True, errors='ignore')

# Encode Status (for plotting later, not for PCA input)
label_encoder = LabelEncoder()
dataset['Status_Label'] = label_encoder.fit_transform(dataset['Status'])

# Remove 'Status' from features used for PCA
dataset.drop(columns=['Status'], inplace=True)

# Convert everything to numeric (force conversion, coerce errors to NaN)
dataset = dataset.apply(pd.to_numeric, errors='coerce')

# Fill remaining NaNs (you can change this logic if needed)
dataset = dataset.fillna(0)


X = dataset.drop(columns=['Status_Label'])
y = dataset['Status_Label']

selector = VarianceThreshold(threshold=0.5)
X_var = selector.fit_transform(X)
print(f"Reduced features from {X.shape[1]} to {X_var.shape[1]} by variance filtering")

variances = X.var(axis=0)
print(f"Feature variances stats: min={variances.min()}, max={variances.max()}, mean={variances.mean()}")

zero_var_cols = [col for col in X.columns if X[col].var() == 0]
print("Columns with zero variance:", zero_var_cols)
X = X.drop(columns=zero_var_cols)



# Standardize
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_var)

# PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

# Prepare for plotting
pca_df = pd.DataFrame(data=X_pca, columns=['PC1', 'PC2'])
pca_df['Status'] = label_encoder.inverse_transform(y)

# Plot PCA
plt.figure(figsize=(8, 6))
sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='Status', palette='Set1')
plt.title('PCA of Clinical Dataset')
plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)')
plt.legend(title='Status')
plt.grid(True)
plt.tight_layout()
plt.show()
plt.plot(np.cumsum(pca.explained_variance_ratio_))
plt.xlabel('Number of PCs')
plt.ylabel('Cumulative Explained Variance')
plt.title('Scree Plot')
plt.grid()
plt.show()

