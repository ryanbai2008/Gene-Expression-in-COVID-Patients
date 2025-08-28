# `mod_member` and `gene_sig` are your module membership and gene significance values
import matplotlib.pyplot as plt
import pandas as pd
# Example data for plotting
mod_member = pd.read_csv("module_membership.csv")
gene_sig = pd.read_csv("gene_significance.csv")

plt.figure(figsize=(8, 6))
plt.scatter(mod_member, gene_sig)
plt.xlabel("Module Membership (MM)")
plt.ylabel("Gene Significance (GS)")
plt.title("Module Membership vs Gene Significance")
plt.show()
