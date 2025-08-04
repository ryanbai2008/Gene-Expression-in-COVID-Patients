import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#PERFORM DEA FIRST
# Load DESeq2 results (must contain baseMean and log2FoldChange)
df = pd.read_csv("./data2/DESeq2_results.csv")

# Drop rows with missing values
df = df.dropna(subset=["log2FoldChange", "baseMean"])

# Compute log10 of baseMean for A-axis
df['log10(baseMean)'] = np.log10(df['baseMean'] + 1)  # Add 1 to avoid log(0)

# Determine significance (optional: adjust thresholds)
df['significant'] = (df['padj'] < 0.05) & (abs(df['log2FoldChange']) > 1)

# Plot
plt.figure(figsize=(10, 6))
sns.scatterplot(
    data=df,
    x='log10(baseMean)',
    y='log2FoldChange',
    hue='significant',
    palette={True: 'red', False: 'grey'},
    alpha=0.7,
    edgecolor=None
)

plt.axhline(0, linestyle='--', color='black', linewidth=1)
plt.title('MA Plot of Differential Expression')
plt.xlabel('log₁₀(Base Mean Expression)')
plt.ylabel('log₂(Fold Change)')
plt.legend(title='Significant', loc='upper right')
plt.tight_layout()
plt.show()
