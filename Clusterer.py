import math
import numpy as np
from sklearn.cluster import KMeans
import Domain
import DomainBuilder as dom_build

class Clusterer:
    def __init__(self, max_k, params, rotation_vectors, residues_1, residues_2):
        self.max_k: int = max_k
        self.params: {} = params
        self.rotation_vectors = rotation_vectors
        self.residues_1 = residues_1
        self.residues_2 = residues_2
        self.atom_coordinates_1 = []
        self.atom_coordinates_2 = []
        self.k_means_results = None
        self.segments = {}
        self.residue_atoms = []
        self.domains_1 = []
        self.domains_2 = []

    def cluster(self):
        num_iters = 50
        # num_rotation_vecs = self.rotation_vectors.shape[0]
        # num_rotation_vec_axis = self.rotation_vectors.shape[1]
        current_k = 3
        # finished = False
        # while (not finished) and current_k < self.max_k:
        # while current_k < self.max_k:
        while current_k < 4:
            print(f"current_k = {current_k}")
            # KMeans the rotation vectors to obtain k number of clusters
            self.k_means_results = self.calc_k_means(current_k, num_iters)
            # self.print_labels()
            # Check if segments are
            # if not self.check_cluster_size():
            #     print("A cluster contains too few residues")
            #     break
            # Obtain the segments from the KMeans results
            self.segments = self.determine_segments(current_k)
            # self.print_segments()
            print("Protein 1")
            self.domains_1 = dom_build.build_domains(self.residues_1, self.segments, int(self.params["domain"]))
            # print(self.domains_1)
            print("Protein 2")
            self.domains_2 = dom_build.build_domains(self.residues_2, self.segments, int(self.params["domain"]))
            # print(self.domains_2)
            current_k += 1

    def calc_k_means(self, k, iters):
        k_means = KMeans(n_clusters=k, random_state=0, n_init="auto", max_iter=iters).fit(self.rotation_vectors)
        return k_means

    def determine_segments(self, k):
        """
        Finds segments in the array of K-Means cluster IDs where the segments contain the same cluster IDs in a row.
        :param k:   Number of clusters
        :return:    A dictionary where the keys are the cluster IDs and the values are a list of tuples (start, end)
                    where start is the index in k_mean_labels where the segment starts and end is the index in
                    k_mean_labels where the segment ends.
        """
        k_mean_labels = self.k_means_results.labels_
        current_element_to_check = k_mean_labels[0]
        start_index = 0
        end_index = 0
        # This line of code does not work as it will end up appending to all lists
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

    # def check_cluster_size(self):
        # unique, counts = np.unique(self.k_means_results.labels_, return_counts=True)
        # hits = np.where(counts < int(self.params["domain"]))
        # if len(hits[0]) > 1:
        #     return False
        # return True

        # count = 0
        # for cluster_segments in self.segments.values():
        #     cluster_size_reached = False
        #     for segment in cluster_segments:
        #         count += segment[1] + 1 - segment[0]
        #         if count >= int(self.params["domain"]):
        #             count = 0
        #             cluster_size_reached = True
        #             break
        #     if not cluster_size_reached:
        #         return False
        # return True

    def remove_tiny_clusters(self):
        unique, counts = np.unique(self.k_means_results.labels_, return_counts=True)
        print(f"Unique = {unique}")
        print(f"Counts = {counts}")

    def print_labels(self):
        print("Printing Labels...")
        print(f"Length = {len(self.k_means_results.labels_)}")
        print(f"Labels = {[self.k_means_results.labels_[i] for i in range(len(self.k_means_results.labels_))]}")

    def print_segments(self):
        print("Printing segments...")
        count = 0
        for k, v in self.segments.items():
            print(f"Cluster {k}")
            print(f"Values {v}")
            for i in v:
                print(self.k_means_results.labels_[i[0]:i[1]+1])
                count += len(self.k_means_results.labels_[i[0]:i[1]+1])
        print(f"total = {count}")
        return

    # def get_unit_vectors(self):
    #     for i in range(self.rotation_vectors.shape[0]):
    #         vec = self.rotation_vectors[i]
    #         norm_vec = vec / np.linalg.norm(vec)
    #         self.unit_vectors[i] = norm_vec

