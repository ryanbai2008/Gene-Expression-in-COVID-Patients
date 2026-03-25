import pandas as pd
from scipy.stats import chi2_contingency

# loads data
eigengenes = pd.read_csv("module_eigengenes.csv", index_col=0) #change links of these
metadata = pd.read_csv("traits.csv")

# assigns each sample to a cluster
eigengenes_numeric = eigengenes.select_dtypes(include=['float64', 'int64'])
cluster_assignments = eigengenes_numeric.idxmax(axis=1)  # pick module with highest eigengene value
cluster_assignments = pd.DataFrame({
    "Sample": eigengenes_numeric.index,
    "Cluster": cluster_assignments
})

# prepare metadata
metadata = metadata.rename(columns={
    metadata.columns[0]: "Sample",   # patient/sample ID column
    metadata.columns[1]: "Status"    # infection status column
})

# merge clusters + metadata 
merged = pd.merge(cluster_assignments, metadata, on="Sample", how="inner")

# creates contingency table 
contingency_table = pd.crosstab(merged["Cluster"], merged["Status"])

# performs Chi-Square test
chi2, p_value, dof, expected = chi2_contingency(contingency_table)

# prints results
print("\n Contingency Table (Cluster vs Infection Status)")
print(contingency_table)

print("\nResults")
print(f"Chi-Square Statistic: {chi2:.4f}")
print(f"Degrees of Freedom: {dof}")
print(f"P-Value: {p_value:.6f}")

if p_value < 0.05:
    print("\n Clusters are GOOD; strong association between clusters and infection status.")
else:
    print("\n Clusters are BAD; no significant association between clusters and infection status.")
