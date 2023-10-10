import numpy as np
import gemmi
from sklearn.cluster import KMeans
import Protein
import Domain
import DomainBuilder as dom_build


class Clusterer:
    def __init__(self, max_k: int, domain_size: int, ratio: float, rotation_vectors: np.array, protein_1: Protein, protein_2: Protein):
        self.max_k: int = max_k
        self.min_domain_size: int = domain_size
        self.min_ratio: float = ratio
        self.rotation_vectors: np.array = rotation_vectors
        self.protein_1: Protein = protein_1
        self.protein_2: Protein = protein_2
        self.k_means_results = None
        self.segments = {}
        self.domains = []
        self.fixed_domain = None
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
            temp_segments = self.determine_segments(current_k)
            # self.print_segments()
            print("Protein 1")
            temp_domains = dom_build.build_domains(self.protein_1.slide_window_residues, self.segments, self.min_domain_size)
            print(temp_domains)
            # print("Protein 2")
            # temp_2_domains = dom_build.build_domains(self.residues_2, self.segments, self.min_domain_size)
            # print(temp_2_domains)
            if not (self.check_domain_sizes(temp_domains)):
                break
            temp_fixed_domain = self.domain_connectivity(temp_domains)
            protein_fit = self.whole_protein_domain_best_fit(temp_domains)
            # self.calc_ext_int_ratio(temp_1_domains, self.residues_1, self.residues_2, temp_fixed_domain)
            # self.protein_1_domains = temp_1_domains
            # self.protein_2_domains = temp_2_domains
            self.segments = temp_segments
            self.domains = temp_domains
            self.fixed_domain = temp_fixed_domain
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

    def domain_connectivity(self, domains):
        """
        Takes the domains of one of the protein conformation and determines which domain has the most number of domains
        connected to it.
        :return:
        """
        num_domains = len(domains)
        if num_domains <= 2:
            chosen_domain = np.argmax([sum(d.segments[:, 1] + 1 - d.segments[:, 0]) for d in domains])
            return chosen_domain
        segments = [np.sort(d.segments, axis=0) for d in domains]
        # print(f"Segments = {segments}")
        # print(f"Before = {domains[0].segments}")
        # print(f"After = {np.sort(domains[0].segments, axis=0)}")
        connectivity = {d.id: np.array([]) for d in domains}
        for curr_d in domains:
            prev_indices = segments[curr_d.id][:, 0] - 1
            next_indices = segments[curr_d.id][:, 1] + 1
            for d in domains:
                if curr_d.id == d.id:
                    continue
                prev_hits = np.in1d(prev_indices, segments[d.id])
                next_hits = np.in1d(next_indices, segments[d.id])
                if np.any([prev_hits, next_hits]):
                    temp = np.append(connectivity[curr_d.id], d.id)
                    connectivity[curr_d.id] = temp
                    temp = np.append(connectivity[d.id], curr_d.id)
                    connectivity[d.id] = temp
        # Can there be multiple max connectivities?
        # https://datagy.io/python-get-dictionary-key-with-max-value/
        connectivity = {key: np.unique(value).size for key, value in connectivity.items()}
        chosen_domain = max(connectivity, key=connectivity.get)
        return chosen_domain

    def whole_protein_domain_best_fit(self, domains: list, residues_1: list, residues_2: list):

        fit = 0.0
        for i in range(len(domains)):
            domain: Domain = domains[i]
            mass = 1 / domain.num_residues
            
        return fit

    def superimpose_chains(self):
        """
        Superimposes each residue of the 2 Protein structures using backbone atoms.
        :return: List of gemmi.SupResult objects containing superimposition information
        """
        # Get the backbone atoms
        backbone_1 = self.protein_1.chain_atoms
        backbone_2 = self.protein_2.chain_atoms
        self.chain_superimpose_result = []
        try:
            # For each residue in the 2 proteins, superimpose the backbone atoms
            for i in range(backbone_1.shape[0]):
                # Get the x, y, and z coordinates of the backbone atoms (N, CA, C) of the specific residue
                pos_1 = [a.pos for a in backbone_1[i][:]]
                pos_2 = [a.pos for a in backbone_2[i][:]]
                # Superimpose and append
                self.chain_superimpose_result.append(gemmi.superpose_positions(pos_1, pos_2))
        except Exception as e:
            print(e)

    def calc_ext_int_ratio(self, domains: Domain, backbone_chain_1, backbone_chain_2, fixed_domain_id):
        """
        Calculates the ratio of the inter (external) and intra (internal) motions of domain pairs
        :return:
        """
        ratios = []
        for d in domains:
            if d.id == fixed_domain_id:
                continue
            mass = 1 / d.num_residues
            ext = self.calc_external_motion(d, backbone_chain_1, backbone_chain_2, mass)
            int = self.calc_internal_motion(d, backbone_chain_1, backbone_chain_2, mass)
            # ratio = math.sqrt()
            # ratios.append(ratio)
        return

    def calc_external_motion(self, domain, backbone_chain_1, backbone_chain_2, mass):
        motion = 0.0
        coordinates = []
        for s in domain.segments:
            coordinates.extend(backbone_chain_1[s[0]:s[1]])
        return motion

    def calc_internal_motion(self, domain, backbone_chain_1, backbone_chain_2, mass):
        motion = 0.0
        return motion

    def check_domain_sizes(self, domains):
        return False if any(d.num_residues < self.min_domain_size for d in domains) else True

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

