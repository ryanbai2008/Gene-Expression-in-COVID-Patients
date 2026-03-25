import pandas as pd 
import json

data = pd.read_csv(r"data2\normalized_counts_DESeq2(2).csv")
with open(r"data\Genes used for analysis (TOM).json", 'r') as file:
    desired_gene_names = json.load(file)
all_gene_names = data.columns

# simple filter
filtered = data[data['Unnamed: 0'].isin(desired_gene_names)].reset_index(drop=True)
# get the set of gene IDs present
present = set(data['Unnamed: 0'].astype(str))  # cast to str to be robust

# find missing genes
missing = [g for g in desired_gene_names if str(g) not in present]
if missing:
    print(f"{len(missing)} desired genes not found. Example missing: {missing[:10]}")

# keep only those desired genes that are present, in the same order
keep_in_order = [g for g in desired_gene_names if str(g) in present]

# create filtered df in that order
# set index to gene ID, reindex by keep_in_order, then reset index
filtered_ordered = (data.set_index('Unnamed: 0', drop=False)
                        .loc[keep_in_order]
                        .reset_index(drop=True))
# or, if you used filtered_ordered:
filtered_ordered.to_csv(r"data\filtered_genes_ordered.csv", index=False)

# data format after data.head(3)
#         Unnamed: 0       HC22      HC41      HC43      HC44      HC47  ...   nCIDN001   nCIDN014    nCNE106    nCNE109    nCNE113   nCNE118
# 0  ENSG00000000419   8.401109  8.359582  8.483691  8.452010  8.327780  ...   8.884206   8.872143   8.298127   8.521417   8.428349  8.214371 
# 1  ENSG00000000457  10.039614  9.819373  9.973598  9.904801  9.861249  ...  10.273952  10.510449  10.067713  10.104738  10.497442  9.776683 

# [3 rows x 251 columns]