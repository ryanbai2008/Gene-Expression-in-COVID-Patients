import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# Step 1: Load the edge CSV file (it contains the precomputed edges and weights)
edges_df = pd.read_csv("module_network_edges_.csv")  # Adjust path to your edge file

# Step 2: Create a graph
G = nx.Graph()

# Step 3: Add edges from the edge CSV file (assumes columns: Gene1, Gene2, Weight)
for _, row in edges_df.iterrows():
    G.add_edge(row['Gene1'], row['Gene2'], weight=row['Weight'])

# Step 4: Draw the network (can adjust parameters for better layout)
plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G, k=0.15, iterations=20)
nx.draw(G, pos, with_labels=True, node_size=100, font_size=8, node_color='skyblue', edge_color='gray')
plt.title("Gene Network from Edge CSV")
plt.show()
