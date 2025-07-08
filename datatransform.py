import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import re
import os
os.environ['R_HOME'] = r'C:\Program Files\R\R-4.5.1'

os.environ['PATH'] += r';C:\Program Files\R\R-4.5.1\bin\x64'

import rpy2.robjects as ro
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

pandas2ri.activate()


# Load DESeq2 package in R
robjects.r('''
    library(DESeq2)
''')

df = pd.read_csv("data/GSE217948_combined_expression_matrix(1).csv", index_col=0)

df = df.apply(pd.to_numeric, errors='coerce')

# Drop genes with all NaNs (non-expressed or invalid rows)
df.dropna(axis=0, how='all', inplace=True)

min_counts = 20
min_samples = int(0.5 * df.shape[1])  # 80% of samples

# Filter genes
filtered_df = df[df.gt(min_counts).sum(axis=1) >= min_samples]


# Check total expression per sample
sample_sums = filtered_df.sum(axis=0)

# Optional: visualize total counts to see if outliers are obvious
plt.figure(figsize=(10, 4))
sns.boxplot(data=sample_sums)
plt.title("Total Gene Expression per Sample")
plt.ylabel("Total Counts")
plt.xticks([])
plt.show()

# Remove samples with extreme expression (e.g., beyond 3 std dev)
z_scores = (sample_sums - sample_sums.mean()) / sample_sums.std()
non_outliers = z_scores.abs() < 3
filtered_df = filtered_df.loc[:, non_outliers]
print("Outlier samples:")
print(sample_sums[~non_outliers])


print(f"Filtered from {df.shape[0]} genes to {filtered_df.shape[0]} genes")

def run_deseq2(counts_df, meta_df):
    from rpy2.robjects import r, pandas2ri
    pandas2ri.activate()
    r_counts = pandas2ri.py2rpy(counts_df)
    r_meta = pandas2ri.py2rpy(meta_df)
    
    robjects.globalenv['countData'] = r_counts
    robjects.globalenv['colData'] = r_meta
    
    robjects.r('''
        dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
        dds <- DESeq(dds)
        normalized_counts <- counts(dds, normalized=TRUE)
    ''')
    
    normalized = robjects.r['normalized_counts']
    return pandas2ri.rpy2py(normalized)

def normalize_sample_name(name):
    match = re.match(r"^(n\w+?)(d\d+)?$", name)
    if match:
        return match.group(1)
    return name

metadata = pd.read_excel("data/data.xlsx")
metadata = metadata.set_index('Patient_ID')  # Ensure 'Patient_ID' is the index
print(metadata['Status'].unique())
filtered_df.columns = filtered_df.columns.to_series().apply(normalize_sample_name)


norm_counts = run_deseq2(filtered_df, metadata)
