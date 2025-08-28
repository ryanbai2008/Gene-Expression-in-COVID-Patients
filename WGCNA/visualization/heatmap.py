import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# loading data
df = pd.read_csv("./data2/normalized_counts_DESeq2(2).csv")
df = df.set_index(df.columns[0])
df['variance'] = df.var(axis=1)
df_top = df.sort_values(by='variance', ascending=False).head(50)
df_top = df_top.drop(columns=['variance'])
sns.set(style="whitegrid")
sns.clustermap(
 df_top,
 cmap="viridis",
 figsize=(12, 10),
 standard_scale=0 # 0 = normalize genes (rows)
)
plt.title("Heatmap of Top 50 Most Variable Genes")
plt.show()