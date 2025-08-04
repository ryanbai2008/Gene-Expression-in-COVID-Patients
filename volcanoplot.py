import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Load DESeq2 results (PERFORM DEA FIRST, CHANGE FILE)
df = pd.read_csv("./data2/DESeq2_results.csv")  # Make sure this file has log2FoldChange and padj

# Drop rows with NA p-values
df = df.dropna(subset=["log2FoldChange", "padj"])

# Add a column for significance
df['-log10(padj)'] = -np.log10(df['padj'])
df['significant'] = (df['padj'] < 0.05) & (abs(df['log2FoldChange']) > 1)

# Plot
plt.figure(figsize=(10, 8))
sns.scatterplot(
    data=df,
    x='log2FoldChange',
    y='-log10(padj)',
    hue='significant',
    palette={True: 'red', False: 'grey'},
    alpha=0.7,
    edgecolor=None
)

plt.axhline(-np.log10(0.05), linestyle='--', color='black', linewidth=1)
plt.axvline(-1, linestyle='--', color='black', linewidth=1)
plt.axvline(1, linestyle='--', color='black', linewidth=1)

plt.title('Volcano Plot of Differentially Expressed Genes')
plt.xlabel('Log₂ Fold Change')
plt.ylabel('-Log₁₀ Adjusted p-value')
plt.legend(title='Significant', loc='upper left')
plt.tight_layout()
plt.show()
