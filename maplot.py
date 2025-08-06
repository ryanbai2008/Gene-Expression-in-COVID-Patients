import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load DE results
df = pd.read_csv("./data2/limma_all_results.csv")  

# Drop rows with missing values
df = df.dropna(subset=["logFC", "AveExpr", "adj.P.Val"])

# Add significance flag
df['significant'] = (df['adj.P.Val'] < 0.05) & (abs(df['logFC']) > 0.58)

# Plot
plt.figure(figsize=(10, 6))
plt.scatter(df['AveExpr'], df['logFC'], c=df['significant'].map({True: 'red', False: 'grey'}), alpha=0.6, edgecolors='none')

plt.axhline(0, color='black', linestyle='--', linewidth=1)

plt.xlabel('Average Log₂ Expression (A)')
plt.ylabel('Log₂ Fold Change (M)')
plt.title('MA Plot (Limma)')

# Optional: show number of significant genes
n_sig = df['significant'].sum()
plt.text(df['AveExpr'].min(), df['logFC'].max(), f'Significant DEGs: {n_sig}', color='red')

plt.tight_layout()
plt.show()

print("Total genes:", len(df))
print("Significant DEGs:", df['significant'].sum())
