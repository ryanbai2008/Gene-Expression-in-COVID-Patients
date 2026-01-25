import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

os.environ["OMP_NUM_THREADS"] = "2"

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from data.getData import *

from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

# data!!
X, y = getData()

# scaling
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

#range of k values
k_vals = range(1,11)
wcss = []

#find wcss
for k in k_vals:
    kmeans = KMeans(n_clusters=k, random_state=42, n_init='auto')
    kmeans.fit(X_scaled)
    wcss.append(kmeans.inertia_)

# getting elbow
def findElbow(wcss):
    k_vals = np.arange(1,len(wcss)+1)
    lineStart = np.array([k_vals[0], wcss[0]])
    lineEnd = np.array([k_vals[-1], wcss[-1]])

    dist = []
    for i in range(len(k_vals)):
        pt = np.array([k_vals[i], wcss[i]])
        #dist from pt to line
        num = np.abs(np.cross(lineEnd-lineStart, lineStart-pt))
        den = np.linalg.norm(lineEnd-lineStart)
        dist.append(num/den)
    
    elbowK = k_vals[np.argmax(dist)]
    return elbowK

elbow_k = findElbow(wcss)

#plot wcss
plt.figure(figsize=(8,5))
plt.plot(k_vals, wcss, marker='o')
plt.scatter(elbow_k, wcss[elbow_k-1], color='green')
plt.title('elbows')
plt.xlabel('cluster nums')
plt.ylabel('wcss')
plt.grid(True)
plt.tight_layout()
plt.show()
