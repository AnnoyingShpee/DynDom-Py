import math

import numpy as np
# from Collections import Counter
import sklearn.cluster._kmeans
from sklearn.cluster import KMeans


class Clusterer:
    def __init__(self, max_k, params, rotation_vectors):
        self.max_k: int = max_k
        self.params: {} = params
        self.rotation_vectors = rotation_vectors
        self.segments = {}
        self.segment_mid_points = {}

        # self.unit_vectors = np.empty(shape=rotation_vectors.shape)

    def cluster(self):
        num_iters = 50
        # num_rotation_vecs = self.rotation_vectors.shape[0]
        # num_rotation_vec_axis = self.rotation_vectors.shape[1]
        current_k = 1
        finished = False
        while not finished:
            results: sklearn.cluster._kmeans.KMeans = self.calc_k_means_sklearn(current_k, num_iters)
            segments = self.determine_segments(results.labels_, current_k)
            # self.print_segments(results.labels_, segments)
            segment_mid_points = self.get_segment_middle_element(segments)
            if current_k >= 3:
                finished = True
                self.segments = segments
                self.segment_mid_points = segment_mid_points
            else:
                current_k += 1

    def determine_segments(self, k_mean_labels, k):
        """
        Finds segments in the array of K-Means cluster IDs where the segments contain the same cluster IDs in a row.
        :param k_mean_labels: An array of K-Means cluster IDs
        :param k:   Number of clusters
        :return:    A dictionary where the keys are the cluster IDs and the values are a list of tuples (start, end)
                    where start is the index in k_mean_labels where the segment starts and end is the index in
                    k_mean_labels where the segment ends.
        """
        current_element_to_check = k_mean_labels[0]
        start_index = 0
        end_index = 0
        # This does not work as it will end up appending to all lists
        # segment_indices = dict.fromkeys(range(0, k), [])
        # This is the correct way to initialise a dictionary of lists
        segment_indices = {key: [] for key in range(k)}
        for i in range(len(k_mean_labels)):
            if k_mean_labels[i] != current_element_to_check:
                segment_indices[current_element_to_check].append((start_index, end_index))
                current_element_to_check = k_mean_labels[i]
                start_index = i
            if i == len(k_mean_labels) - 1:
                segment_indices[current_element_to_check].append((start_index, i))
            end_index = i
        return segment_indices

    def get_segment_middle_element(self, segments):
        """
        Gets the middle element of the segments
        :param segments:
        :return: Dictionary of lists of indexes of k_mean_labels
        """
        segment_mid_index = {key: [] for key in range(len(segments))}
        for key, value in segments.items():
            for v in value:
                index = (v[0] + v[1]) // 2
                segment_mid_index[key].append(index)
        return segment_mid_index

    def calc_k_means_sklearn(self, k, iters):
        k_means = KMeans(n_clusters=k, random_state=0, n_init="auto", max_iter=iters).fit(self.rotation_vectors)
        # print(f"Cluster labels = {k_means.labels_}")
        # print(f"Cluster centers = {k_means.cluster_centers_}")
        return k_means

    # def check_criterion(self, k):
    #     return True

    def print_segments(self, k_mean_labels, segments):
        print(f"labels = {[k_mean_labels[i] for i in range(len(k_mean_labels))]}")
        print(f"Length = {len(k_mean_labels)}")
        count = 0
        for k, v in segments.items():
            print(f"Cluster {k}")
            print(f"Values {v}")
            for i in v:
                print(k_mean_labels[i[0]:i[1]+1])
                count += len(k_mean_labels[i[0]:i[1]+1])
        print(f"total = {count}")
        return

    # def get_unit_vectors(self):
    #     for i in range(self.rotation_vectors.shape[0]):
    #         vec = self.rotation_vectors[i]
    #         norm_vec = vec / np.linalg.norm(vec)
    #         self.unit_vectors[i] = norm_vec