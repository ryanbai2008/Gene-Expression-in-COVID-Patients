import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd

# Step 1: Load the TOM dissimilarity matrix (adjust the path if necessary)
TOM_diss_df = pd.read_csv("./data/tom_matrix_values.csv", index_col=0)

# Convert dissimilarity to similarity
TOM_sim = 1 - TOM_diss_df.values
genes = TOM_diss_df.columns

# Step 2: Load the hub genes from your hub_genes.csv file
hub_genes_df = pd.read_csv("./data/hub_genes_0.5.csv")
hub_genes = hub_genes_df['Hub Genes'].tolist()  # Extract the hub genes list

# Step 3: Filter TOM similarity matrix to include only hub genes
# Find the indices of the hub genes in the TOM similarity matrix
hub_genes_indices = [genes.get_loc(gene) for gene in hub_genes if gene in genes]

# Create a submatrix for the hub genes
TOM_submatrix = TOM_sim[hub_genes_indices, :][:, hub_genes_indices]

# Step 4: Create a graph for the hub gene subnetwork using the submatrix
G_hub = nx.Graph()

# Add nodes for each hub gene
G_hub.add_nodes_from(hub_genes)

# Add edges based on the TOM similarity (thresholding)
threshold = 0.99995  # Adjust this threshold based on your needs
for i in range(len(hub_genes)):
    for j in range(i + 1, len(hub_genes)):
        if TOM_submatrix[i, j] > threshold:
            G_hub.add_edge(hub_genes[i], hub_genes[j], weight=TOM_submatrix[i, j])

# Step 5: Draw the subnetwork of hub genes
plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G_hub, k=0.2, iterations=50)  # Adjust layout parameters for better spacing
nx.draw(G_hub, pos, with_labels=True, node_size=100, font_size=8, node_color='skyblue', edge_color='gray')
plt.title("Hub Gene Subnetwork (Topological Overlap)")
plt.show()
