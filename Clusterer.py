import sklearn.cluster._kmeans
from sklearn.cluster import KMeans



def calc_k_means_sklearn(data, k):
    k_means: sklearn.cluster._kmeans.KMeans = KMeans(n_clusters=k, random_state=0, n_init="auto").fit(data)
    print(f"Cluster labels = {k_means.labels_}")
    print(f"Cluster centers = {k_means.cluster_centers_}")
    return k_means.labels_


