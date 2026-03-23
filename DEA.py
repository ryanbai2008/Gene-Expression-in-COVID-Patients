import pandas as pd
import numpy as np
import re
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.feature_selection import VarianceThreshold


# Load expression matrix (genes as rows, samples as columns)
df = pd.read_csv("./data2/normalized_counts_DESeq2(2).csv", index_col=0)

# Transpose: rows = samples, columns = genes
X = df.T.copy()

# Normalize sample names (remove trailing d#)
def normalize_sample_name(name):
    return re.sub(r'd\d+$', '', name)

X['Normalized'] = X.index.to_series().apply(normalize_sample_name)

# Load metadata, set Patient_ID as index
metadata = pd.read_excel("data/data.xlsx")
metadata = metadata.set_index('Patient_ID')

# Merge metadata into expression data by normalized sample names
X = X.merge(metadata, left_on='Normalized', right_index=True, how='inner')

# Extract group labels
group_labels = X['Status']

# Drop helper columns
X.drop(columns=['Normalized', 'Status'], inplace=True)

# Convert expression data to numeric and fill NaNs with 0
X = X.apply(pd.to_numeric, errors='coerce').fillna(0)



# Filter for two groups only: e.g., "COVID-19" vs "Control"
group1 = "COVID-19"
group2 = "Control"
valid_samples = group_labels[(group_labels == group1) | (group_labels == group2)].index
X_filtered = X.loc[valid_samples]
labels_filtered = group_labels.loc[valid_samples]


selector = VarianceThreshold(threshold=0.37)
X_highvar = selector.fit_transform(X_filtered)
selected_genes = X_filtered.columns[selector.get_support()]
X_filtered = pd.DataFrame(X_highvar, index=X_filtered.index, columns=selected_genes)

# Run linear model gene by gene
results = []
for gene in X_filtered.columns:
    y = X_filtered[gene]
    # Design matrix with intercept and binary group indicator
    X_design = pd.DataFrame({
        'intercept': 1,
        'group': labels_filtered.map({group1: 1, group2: 0})
    })
    model = sm.OLS(y, X_design).fit()
    log2fc = model.params['group']
    pval = model.pvalues['group']
    results.append({'gene': gene, 'log2FoldChange': log2fc, 'pval': pval})

dea_df = pd.DataFrame(results)

print("Checking p-values before correction...")
print(dea_df['pval'].isna().sum(), "NaN p-values found out of", len(dea_df))
print(dea_df['pval'].head())


# Drop rows with NaN p-values before multiple testing correction
dea_df_clean = dea_df.dropna(subset=['pval']).copy()

print("Running multiple testing correction on", len(dea_df_clean), "valid p-values")


padj_corrected = multipletests(dea_df_clean['pval'], method='fdr_bh')[1]
dea_df_clean['padj'] = padj_corrected

# Merge back corrected p-values to original table
dea_df = dea_df.drop(columns=['padj'], errors='ignore')  # In case padj exists already
dea_df = dea_df.merge(dea_df_clean[['gene', 'padj']], on='gene', how='left')

# Sort by adjusted p-value
dea_df = dea_df.sort_values('padj')

# Volcano plot
dea_df['-log10(padj)'] = -np.log10(dea_df['padj'])

plt.figure(figsize=(10,6))
sns.scatterplot(data=dea_df, x='log2FoldChange', y='-log10(padj)',
                hue=dea_df['padj'] < 0.05,
                palette={True: 'red', False: 'gray'},
                legend=False)
plt.axhline(-np.log10(0.05), linestyle='--', color='black')
plt.axvline(1, linestyle='--', color='black')
plt.axvline(-1, linestyle='--', color='black')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-log10 Adjusted P-value')
plt.title('Volcano Plot of Differential Expression')
plt.tight_layout()
plt.show()

# Export results
dea_df.to_csv("dea_results.csv", index=False)

print(dea_df[['log2FoldChange', 'padj']].describe())
print(f"Number of significant genes (padj < 0.05): {(dea_df['padj'] < 0.05).sum()}")

print(dea_df[(abs(dea_df['log2FoldChange']) > 1) & (dea_df['padj'] < 0.05)])
print(f"Number of NaN p-values: {dea_df['pval'].isna().sum()}")

#print(dea_df.head(10))
