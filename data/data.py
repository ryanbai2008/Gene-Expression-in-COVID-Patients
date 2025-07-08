import pandas as pd
import os
from glob import glob

# Path to extracted files
data_folder = './data/GSE217948_RAW'

# Match all text files
all_files = glob(os.path.join(data_folder, "*.txt"))

expression_matrix = None

for file in all_files:
    # Extract sample name from filename
    sample_name = os.path.basename(file).split('_')[1].replace('.counts', '').replace('.txt', '')

    # Read file
    df = pd.read_csv(file, sep="\t", header=None)
    df.columns = ['Gene', sample_name]

    # Merge
    if expression_matrix is None:
        expression_matrix = df
    else:
        expression_matrix = pd.merge(expression_matrix, df, on='Gene', how='outer')

# Set 'Gene' as index
expression_matrix.set_index('Gene', inplace=True)

# Fill missing values with 0
expression_matrix = expression_matrix.fillna(0)

# Save to CSV
expression_matrix.to_csv("GSE217948_combined_expression_matrix(1).csv")

print("✅ Expression matrix created with shape:", expression_matrix.shape)
print(expression_matrix)
