import pandas as pd 

KMeans_df = pd.read_csv("KMeans_clusters.csv")
kernel_KMeans_df = pd.read_csv("kernel_KMeans_clusters.csv")
data_df = pd.read_excel("data/data.xlxs")
clusters = [KMeans_df, kernel_KMeans_df]

for cluster in clusters:
  df = data.df.merge(KMeans_df, on="Patient")
  from sklearn.metrics import silhouette_score
  features = df[["Age","Weight"]]
  labels = df["ClusterID"]
  score = silhouette_score(features, labels, metric="euclidean")
  print("Silhouette Score:", score)
