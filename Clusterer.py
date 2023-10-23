import math
import numpy as np
import gemmi
from sklearn.cluster import KMeans
from copy import deepcopy
import Protein
import Domain
from scipy.linalg import svd
import DomainBuilder as dom_build


backbone_atoms = ["N", "CA", "C"]


class Clusterer:
    def __init__(self, max_k: int, params: dict, rotation_vectors: np.array, protein_1: Protein, protein_2: Protein):
        self.max_k: int = max_k
        self.min_domain_size: int = int(params["domain"])
        self.min_ratio: float = float(params["ratio"])
        self.atoms_type = params["atoms"]
        self.num_atoms = 3 if params["atoms"] == "backbone" else 1
        self.rotation_vectors: np.array = rotation_vectors
        self.protein_1: Protein = protein_1
        self.protein_2: Protein = protein_2
        self.k_means_results = None
        self.segments = {}
        self.domains = []
        self.fixed_domain = None

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
            # Obtain the segments from the KMeans results. The segments' indices are for the slide windowed residues.
            temp_segments = self.determine_segments(current_k)
            temp_domains = dom_build.build_domains(self.protein_1.slide_window_residues, temp_segments, self.min_domain_size)
            # print(temp_domains)
            # print("Protein 2")
            # temp_2_domains = dom_build.build_domains(self.residues_2, self.segments, self.min_domain_size)
            # print(temp_2_domains)
            if not (self.check_domain_sizes(temp_domains)):
                print("Size Break")
                break
            temp_fixed_domain = self.domain_connectivity(temp_domains)
            # Perform mass-weighted whole-protein best fit to obtain a new set of coordinates for Protein 1
            r, transformed_1_on_2, transformed_2_on_1 = self.whole_protein(temp_domains)
            # self.print_coords(transformed_1_on_2, transformed_2_on_1)
            # print(transformed)
            ints, rs = self.calc_int(temp_domains, transformed_2_on_1)
            exts = self.calc_ext(temp_domains, rs)
            print(exts)
            print(ints)
            ratio = self.calc_ext_int_ratio(exts, ints, temp_domains)
            print(ratio)

            if ratio < self.min_ratio:
                print("Ratio Break")
                break
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
        # print(len(k_mean_labels))
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

    def whole_protein(self, domains):
        """
        Performs a mass-weighted whole-protein best fit.
        :param domains:
        :return r: A SupResult object containing RMSD, Center 1, Center 2, Rotation Matrix, and Translation Vector
        :return chain_1_copy: The slide window residue chain of Protein 1 after transformation to Protein 2's position
        :return chain_2_copy: The slide window residue chain of Protein 2 after transformation to Protein 1's position
        """
        slide_indices = self.protein_1.slide_window_residues_indices
        chain_1_copy = self.protein_1.residue_span
        chain_2_copy = self.protein_2.residue_span
        for i in range(len(chain_1_copy) - slide_indices[1]):
            chain_1_copy.__delitem__(-1)
            chain_2_copy.__delitem__(-1)
        for i in range(slide_indices[0]):
            chain_1_copy.__delitem__(0)
            chain_2_copy.__delitem__(0)
        # Reset the structures after deletion as the ResidueSpans use referencing. This will "force" chain_2_copy to
        # become a separate object that doesn't reference the Protein's ResidueSpan object
        self.protein_1.recreate_structure()
        self.protein_2.recreate_structure()
        weights = []
        coords_1 = []
        coords_2 = []
        for d in range(len(domains)):
            for s in domains[d].segments:
                for si in range(s[0], s[1] + 1):
                    for a in backbone_atoms:
                        coords_1.append(chain_1_copy[si].sole_atom(a).pos)
                        coords_2.append(chain_2_copy[si].sole_atom(a).pos)
                        weights.append(1/domains[d].num_residues)
        r_2: gemmi.SupResult = gemmi.superpose_positions(coords_1, coords_2, weights)
        r_1: gemmi.SupResult = gemmi.superpose_positions(coords_2, coords_1, weights)
        chain_1_copy.transform_pos_and_adp(r_1.transform)
        chain_2_copy.transform_pos_and_adp(r_2.transform)
        # print(f"Count = {r_2.count}")
        # print(f"Center 1 = {r_2.center1}")
        # print(f"Center 2 = {r_2.center2}")
        # print(f"RMSD = {r_2.rmsd}")
        # print(f"Rotation = {r_2.transform.mat}")
        # print(f"Translation = {r_2.transform.vec}")
        return r_2, chain_1_copy, chain_2_copy

    def superimpose(self, domains, residue_span, protein=1):
        """

        :param domains:
        :param residue_span: A ResidueSpan object
        :param protein: The Protein to be used to superimpose onto the residue_span
        :return:
        """
        slide_indices = self.protein_1.slide_window_residues_indices
        results = []
        if protein == 1:
            chain = self.protein_1.residue_span
            for i in range(len(chain) - slide_indices[1]):
                chain.__delitem__(-1)
            for i in range(slide_indices[0]):
                chain.__delitem__(0)
            self.protein_1.recreate_structure()
        else:
            chain = self.protein_2.residue_span
            for i in range(len(chain) - slide_indices[1]):
                chain.__delitem__(-1)
            for i in range(slide_indices[0]):
                chain.__delitem__(0)
            self.protein_2.recreate_structure()
        for d in range(len(domains)):
            bottom_coords = []
            top_coords = []
            for s in domains[d].segments:
                for si in range(s[0], s[1] + 1):
                    for a in backbone_atoms:
                        top_coords.append(chain[si].sole_atom(a).pos)
                        bottom_coords.append(residue_span[si].sole_atom(a).pos)
            r = gemmi.superpose_positions(bottom_coords, top_coords)
            results.append(r)
        return results

    def calc_int(self, domains, transformed_whole_protein):
        """

        :param domains:
        :param transformed_whole_protein: The slide window residue chain of protein 2 after superposition onto Protein 1
        :return:
        """
        int_msf = []
        rs = self.superimpose(domains, transformed_whole_protein, protein=1)
        for r in range(len(rs)):
            int_msf.append((rs[r].rmsd ** 2) * domains[r].num_residues)

        return int_msf, rs

    def calc_ext(self, domains, rs):
        ext_msf = []
        slide_indices = self.protein_1.slide_window_residues_indices
        # rs = self.superimpose(domains, transformed_2, protein=1)
        for r in range(len(rs)):
            original_coords = np.empty(shape=[domains[r].num_residues * self.num_atoms, 3])
            transformed_coords = np.empty(shape=[domains[r].num_residues * self.num_atoms, 3])
            chain = self.protein_1.residue_span
            for i in range(len(chain) - slide_indices[1]):
                chain.__delitem__(-1)
            for i in range(slide_indices[0]):
                chain.__delitem__(0)
            self.protein_1.recreate_structure()
            chain.transform_pos_and_adp(rs[r].transform)
            i = 0
            for s in domains[r].segments:
                for si in range(s[0], s[1] + 1):
                    for a in backbone_atoms:
                        original_coords[i] = np.asarray(self.protein_1.slide_window_residues[si].sole_atom(a).pos.tolist())
                        transformed_coords[i] = np.asarray(chain[si].sole_atom(a).pos.tolist())
                        i += 1
            disp_vecs = (original_coords - transformed_coords) ** 2
            sum_disp = np.sum(np.sum(disp_vecs, axis=1))
            ext_msf.append(sum_disp)
            self.protein_1.recreate_structure()

        return ext_msf

    def calc_ext_int_ratio(self, exts, ints, domains):
        sum_exts = 0
        sum_ints = 0
        for d in range(len(domains)):
            sum_exts += exts[d] / domains[d].num_residues
            sum_ints += ints[d] / domains[d].num_residues
        ratio = math.sqrt(sum_exts/sum_ints)
        return ratio

    def check_domain_sizes(self, domains):
        return False if any(d.num_residues < self.min_domain_size for d in domains) else True

    def print_labels(self):
        print("Printing Labels...")
        print(f"Length = {len(self.k_means_results.labels_)}")
        print(f"Labels = {[self.k_means_results.labels_[i] for i in range(len(self.k_means_results.labels_))]}")

    def print_segments(self, segments):
        print("Printing segments...")
        seg_count = 0
        res_count = 0
        for k, v in segments.items():
            print(f"Cluster {k}")
            print(f"Values {v}")
            seg_count += v.shape[0]
            for i in v:
                # print(self.k_means_results.labels_[i[0]:i[1]+1])
                res_count += len(self.k_means_results.labels_[i[0]:i[1]+1])
        print(f"Seg total = {seg_count}")
        print(f"Res total = {res_count}")
        return

    def print_coords(self, transformed_1, transformed_2):
        for i in range(len(transformed_1)):
            for a in backbone_atoms:
                print(f"==============================================")
                print(f"P1 : {self.protein_1.slide_window_residues[i].sole_atom(a).pos} -> {transformed_1[i].sole_atom(a).pos}")
                print(f"P2 : {self.protein_2.slide_window_residues[i].sole_atom(a).pos} -> {transformed_2[i].sole_atom(a).pos}")
                print(f"==============================================")



