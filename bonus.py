import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import math
from sklearn import datasets
from sklearn.cluster import KMeans
import mykmeanssp

iris = datasets.load_iris()
matrix = iris.data
k_range =list(range(1,11))
inertias = np.zeros(11)

def inertia(k, iris):
    """ Calculate inertia """
    kmeans = KMeans(n_clusters=k, init="k-means++", random_state=0)
    predict = kmeans.fit(iris.data)
    return kmeans.inertia_

for k in k_range:
    inertias[k] = inertia(k, iris)

K = [None, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
plt.plot(K,inertias)
plt.title('Elbow method for selection optimal optimal "K" clusters')
plt.xlabel('K')
plt.ylabel('Average Dispersion')
plt.annotate("Elbow point", (2, inertias[2]), xytext=(3, 200), arrowprops={'arrowstyle': '-|>','color':"k"})
plt.scatter(2, inertias[2], s=500, facecolors='none', edgecolors='k', linestyle="--")
plt.savefig("elbow.png")

