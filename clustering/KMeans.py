# prevent memory leak issues n stuff
import os
os.environ['OMP_NUM_THREADS'] = '1'

import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

# get data, make sure to transpose bc we need it in that format
data_df = pd.read_csv("data2/normalized_counts_DESeq2(2).csv", index_col=False)
data_df = data_df.T
print(data_df.shape)
data_df = data_df.iloc[1:, :]

kmeans = KMeans(n_clusters=5, random_state=10, n_init=3)
kmeans_labels = kmeans.fit_predict(data_df)

print(kmeans_labels.shape)
print(kmeans_labels)
