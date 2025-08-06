import pandas as pd
from mygene import MyGeneInfo

# Load your data
df = pd.read_csv("./data2/limma_significant_DEGs.csv")  # Replace with your actual filename

# Extract the Ensembl IDs (assumes column is called 'Gene')
ensembl_ids = df['gene'].tolist()

# Initialize mygene and query
mg = MyGeneInfo()
query_result = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human')

# Convert query result to DataFrame
mapping_df = pd.DataFrame(query_result)[['query', 'symbol']]
mapping_df.columns = ['Gene', 'GeneSymbol']

# Merge gene symbols into original dataframe
df_merged = df.merge(mapping_df, on='Gene', how='left')

# Save or inspect
df_merged.to_csv("count_matrix_with_symbols.csv", index=False)
