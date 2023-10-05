import math
import numpy as np
from sklearn.cluster import KMeans
import Domain
import DomainBuilder as dom_build


class Clusterer:
    def __init__(self, max_k: int, domain_size: int, rotation_vectors: np.array, residues_1: list, residues_2: list):
        self.max_k: int = max_k
        self.domain_size: int = domain_size
        self.rotation_vectors: np.array = rotation_vectors
        self.residues_1: list = residues_1
        self.residues_2: list = residues_2
        self.atom_coordinates_1 = []
        self.atom_coordinates_2 = []
        self.k_means_results = None
        self.segments = {}
        self.residue_atoms = []
        self.domains_1 = []
        self.domains_2 = []
        self.ext_int_ratio_1: float = 0.0
        self.ext_int_ratio_2: float = 0.0


    def cluster(self):
        num_iters = 50
        current_k = 3
        # finished = False
        # while (not finished) and current_k < self.max_k:
        # while current_k < self.max_k:
        while current_k < 4:
            print(f"current_k = {current_k}")
            # KMeans the rotation vectors to obtain k number of clusters
            self.k_means_results = self.calc_k_means(current_k, num_iters)
            # self.print_labels()
            # Obtain the segments from the KMeans results
            self.segments = self.determine_segments(current_k)
            # self.print_segments()
            print("Protein 1")
            self.domains_1 = dom_build.build_domains(self.residues_1, self.segments, self.domain_size)
            # print(self.domains_1)
            print("Protein 2")
            self.domains_2 = dom_build.build_domains(self.residues_2, self.segments, self.domain_size)
            # print(self.domains_2)
            # if type(self.domains_1) == bool or type(self.domains_2) == bool:
            #     break
            current_k += 1

    def calc_k_means(self, k, iters):
        """
        Performs KMeans on the rotation vectors
        :param k: Number of resulting clusters
        :param iters: Number of maximum iterations
        :return: KMeans results
        """
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
        # Obtains the labels of the KMeans
        k_mean_labels = self.k_means_results.labels_
        # Set the first label as the checker to check the labels
        current_element_to_check = k_mean_labels[0]
        start_index = 0
        end_index = 0
        # This line of code to declare and initialise a dictionary of lists does not work as it will end up
        # appending to all lists.
        # segment_indices = dict.fromkeys(range(0, k), [])
        # This is the correct way to initialise a dictionary of lists
        segment_indices = {key: np.array([[]], dtype=int) for key in range(k)}
        # Iterate through each label
        for i in range(len(k_mean_labels)):
            # When the label is not equal to the checker. It means the segment ends there and the segment's start and
            # end indices are obtained and stored.
            if k_mean_labels[i] != current_element_to_check:
                temp = np.append(segment_indices[current_element_to_check], [[start_index, end_index]],
                                 axis=0 if segment_indices[current_element_to_check].shape[1] > 0 else 1)
                segment_indices[current_element_to_check] = temp
                current_element_to_check = k_mean_labels[i]
                start_index = i
            if i == len(k_mean_labels) - 1:
                temp = np.append(segment_indices[current_element_to_check], [[start_index, end_index]],
                                 axis=0 if segment_indices[current_element_to_check].shape[1] > 0 else 1)
                segment_indices[current_element_to_check] = temp
                # segment_indices[current_element_to_check].append((start_index, i))
            end_index = i
        return segment_indices

    def domain_connectivity(self):
        """
        Takes the domains of one of the protein conformation and determines which domain has the most number of domains
        connected to it.
        :return:
        """
        num_domains = len(self.domains_1)
        chosen_domain = None
        print(self.domains_1)
        # if num_domains <= 2:
        #     chosen_domain = np.argmax([sum(d.segments[:, 1] + 1 - d.segments[:, 0]) for d in domains])
        #     return chosen_domain
        segments = [np.sort(d.segments, axis=0) for d in self.domains_1]
        print(f"Segments = {segments}")
        # print(f"Before = {domains[0].segments}")
        # print(f"After = {np.sort(domains[0].segments, axis=0)}")
        connectivity = {key: np.array([]) for key in range(num_domains)}
        for s in range(num_domains - 1):
            prev_indices = segments[s][:, 0] - 1
            next_indices = segments[s][:, 1] + 1
            for c in range(num_domains):
                if s == c:
                    continue
                prev_hits = np.in1d(prev_indices, segments[c])
                next_hits = np.in1d(next_indices, segments[c])
                if np.any([prev_hits, next_hits]):
                    temp = np.append(connectivity[s], c)
                    connectivity[s] = temp
                    temp = np.append(connectivity[c], s)
                    connectivity[c] = temp
        # Can there be multiple max connectivities?
        # https://datagy.io/python-get-dictionary-key-with-max-value/
        connectivity = {key: np.unique(value).size for key, value in connectivity.items()}
        chosen_domain = max(connectivity, key=connectivity.get)
        return chosen_domain

    def calc_ext_int_ratio(self):
        """

        :return:
        """
        domain_1 = self.domains_1
        domain_2 = self.domains_2
        return

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
        seg_count = 0
        res_count = 0
        for k, v in self.segments.items():
            print(f"Cluster {k}")
            print(f"Values {v}")
            seg_count += v.shape[0]
            for i in v:
                # print(self.k_means_results.labels_[i[0]:i[1]+1])
                res_count += len(self.k_means_results.labels_[i[0]:i[1]+1])
        print(f"Seg total = {seg_count}")
        print(f"Res total = {res_count}")
        return

    # def get_unit_vectors(self):
    #     for i in range(self.rotation_vectors.shape[0]):
    #         vec = self.rotation_vectors[i]
    #         norm_vec = vec / np.linalg.norm(vec)
    #         self.unit_vectors[i] = norm_vec

