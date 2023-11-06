import math
import numpy as np
import gemmi
from sklearn.cluster import KMeans
import Protein
import DomainBuilder as dom_build


class Clusterer:
    def __init__(self, max_k: int, params: dict, rotation_vectors: np.array, protein_1: Protein, protein_2: Protein, main_atoms):
        self.max_k: int = max_k
        self.min_domain_size: int = int(params["domain"])
        self.min_ratio: float = float(params["ratio"])
        self.atoms_type = params["atoms"]
        self.backbone_atoms = main_atoms
        self.num_atoms = len(self.backbone_atoms)
        self.rotation_vectors: np.array = rotation_vectors
        self.protein_1: Protein = protein_1
        self.protein_2: Protein = protein_2
        self.k_means_results = None
        self.segments = {}
        self.domains = []
        # The number determining which domain will be the fixed point
        self.fixed_domain = None
        # List of gemmi.SupResult objects
        self.rs = None

    def cluster(self):
        num_iters = 50
        current_k = 3
        condition_unmet = False
        # while (not finished) and current_k < self.max_k:
        while current_k < self.max_k and not condition_unmet:
            print(f"current_k = {current_k}")
            # KMeans the rotation vectors to obtain k number of clusters
            self.k_means_results = self.calc_k_means(current_k, num_iters)
            # self.print_labels()
            # Obtain the segments from the KMeans results. The segments' indices are for the slide windowed residues.
            temp_segments = self.determine_segments(current_k)
            # self.print_segments(temp_segments)
            if not self.check_cluster_sizes(temp_segments):
                print("Cluster Break")
                condition_unmet = True
            temp_domains_1, cluster_break = dom_build.build_domains(self.protein_1.slide_window_residues, temp_segments, self.min_domain_size)
            # print("Domains 1")
            # self.print_domains(temp_domains_1)
            if cluster_break:
                print("Domain Size Break")
                condition_unmet = True
            temp_fixed_domain_id = self.check_domain_connectivity(temp_domains_1)
            # Perform mass-weighted whole-protein best fit to obtain a new set of Protein 2 coordinates after
            # superposition Protein 2 onto Protein 1
            _, r2, _, transformed_2_on_1 = self.mass_weighted_whole_protein_fit(temp_domains_1)
            # self.print_coords(transformed_1_on_2, transformed_2_on_1)
            ratio_met = self.check_ratios(transformed_2_on_1, temp_domains_1, temp_fixed_domain_id)
            if not ratio_met:
                break
            self.segments = temp_segments
            self.domains = temp_domains_1
            self.fixed_domain = temp_fixed_domain_id
            self.rs = r2
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

    def check_domain_connectivity(self, domains):
        """
        Takes the domains of one of the protein conformation and determines which domain has the most number of domains
        connected to it.
        :return chosen_domain: The ID of the domain with the most number of connected domains
        """
        num_domains = len(domains)
        if num_domains <= 2:
            chosen_domain = np.argmax([sum(d.segments[:, 1] + 1 - d.segments[:, 0]) for d in domains])
            return chosen_domain
        segments = [np.sort(d.segments, axis=0) for d in domains]
        connectivity = {d.domain_id: np.array([]) for d in domains}
        for curr_d in domains:
            prev_indices = segments[curr_d.domain_id][:, 0] - 1
            next_indices = segments[curr_d.domain_id][:, 1] + 1
            for d in domains:
                if curr_d.domain_id == d.domain_id:
                    continue
                prev_hits = np.in1d(prev_indices, segments[d.domain_id])
                next_hits = np.in1d(next_indices, segments[d.domain_id])
                if np.any([prev_hits, next_hits]):
                    temp = np.append(connectivity[curr_d.domain_id], d.domain_id)
                    connectivity[curr_d.domain_id] = temp
                    temp = np.append(connectivity[d.domain_id], curr_d.domain_id)
                    connectivity[d.domain_id] = temp
        # Can there be multiple max connectivities?
        # https://datagy.io/python-get-dictionary-key-with-max-value/
        connectivity = {key: np.unique(value).size for key, value in connectivity.items()}
        print(f"Connectivity {connectivity}")
        chosen_domain = max(connectivity, key=connectivity.get)
        return chosen_domain

    def mass_weighted_whole_protein_fit(self, domains):
        """
        Performs a mass-weighted whole-protein best fit to get a new set of coordinates of a Protein.
        :param domains:
        :return r_1:    A SupResult object containing RMSD, Center 1, Center 2, Rotation Matrix, and Translation Vector
                        of Protein 1
        :return r_2:    A SupResult object containing RMSD, Center 1, Center 2, Rotation Matrix, and Translation Vector
                        of Protein 2
        :return slide_window_1: The slide window residue chain of Protein 1 after fitting (transformation)
                                to Protein 2's position
        :return slide_window_2: The slide window residue chain of Protein 2 after fitting (transformation)
                                to Protein 1's position
        """
        slide_window_1 = self.protein_1.get_slide_window_result()
        slide_window_2 = self.protein_2.get_slide_window_result()
        weights = []
        coords_1 = []
        coords_2 = []
        for d in range(len(domains)):
            for s in domains[d].segments:
                for si in range(s[0], s[1] + 1):
                    for a in self.backbone_atoms:
                        coords_1.append(slide_window_1[si].sole_atom(a).pos)
                        coords_2.append(slide_window_2[si].sole_atom(a).pos)
                        weights.append(1/domains[d].num_residues)

        r_1: gemmi.SupResult = gemmi.superpose_positions(coords_2, coords_1, weights)
        r_2: gemmi.SupResult = gemmi.superpose_positions(coords_1, coords_2, weights)
        slide_window_1.transform_pos_and_adp(r_1.transform)
        slide_window_2.transform_pos_and_adp(r_2.transform)
        # print(f"Count = {r_2.count}")
        # print(f"Center 1 = {r_2.center1}")
        # print(f"Center 2 = {r_2.center2}")
        # print(f"RMSD = {r_2.rmsd}")
        # print(f"Rotation = {r_2.transform.mat}")
        # print(f"Translation = {r_2.transform.vec}")
        return r_1, r_2, slide_window_1, slide_window_2

    def superimpose_domain_to_transformed(self, residue_span, domains, protein=1):
        """
        Superimposes each Protein domain's original residue atoms onto the transformed residue atoms
        :param domains:
        :param residue_span: A ResidueSpan object
        :param protein: The ID of the Protein to be used to superimpose onto residue_span
        :return results: A list of gemmi.SupResult objects
        """
        # List of fitting protein domain chains
        fitting_domain_chains = []
        # List of target protein domain chains
        target_domain_chains = []
        # List of gemmi.SupResult objects
        results = []
        # Get the original chain from the Protein that will be used to superimpose onto the target
        protein: Protein = self.protein_1 if protein == 1 else self.protein_2
        slide_chain = protein.get_slide_window_result()
        for d in range(len(domains)):
            fitting_chain: gemmi.Chain = gemmi.Chain(protein.chain_param)
            target_chain: gemmi.Chain = gemmi.Chain(protein.chain_param)
            # Domain residue atoms as the fitting
            fitting_coords = []
            # The transformed residue atoms as the target
            target_coords = []
            for s in domains[d].segments:
                for si in range(s[0], s[1] + 1):
                    fitting_chain.add_residue(slide_chain[si])
                    target_chain.add_residue(residue_span[si])
                    for a in self.backbone_atoms:
                        fitting_coords.append(slide_chain[si].sole_atom(a).pos)
                        target_coords.append(residue_span[si].sole_atom(a).pos)
            fitting_domain_chains.append(fitting_chain)
            target_domain_chains.append(target_chain)
            # Superpose the domain residue atoms onto the target residue atoms
            r = gemmi.superpose_positions(target_coords, fitting_coords)
            results.append(r)
        return fitting_domain_chains, target_domain_chains, results

    def check_ratios(self, transformed_protein, domains, fixed_domain_id):
        """
        Checks the ratio of internal and external domain movements of each fixed-connected domain pairs.
        :param transformed_protein: The transformed Protein slide window chain fitted onto the other Protein.
        :param domains:
        :param fixed_domain_id:
        :return:
        """
        transformed_domain_chains, protein_domain_chains, results = self.superimpose_domain_to_transformed(transformed_protein, domains)
        fixed_domain_r: gemmi.SupResult = results[fixed_domain_id]
        fixed_domain_num_residues = domains[fixed_domain_id].num_residues
        fixed_domain_int_msf = self.calc_domain_int_msf(fixed_domain_num_residues, fixed_domain_r)
        fixed_domain_ext_msf = self.calc_domain_ext_msf(domains[fixed_domain_id], fixed_domain_r)
        for d in range(len(domains)):
            if d == fixed_domain_id:
                continue
            connected_domain_int_msf = self.calc_domain_int_msf(domains[d].num_residues, results[d])
            connected_domain_ext_msf = self.calc_domain_ext_msf(domains[d], results[d])
            ratio = self.calc_ext_int_ratio(fixed_domain_ext_msf, fixed_domain_int_msf, fixed_domain_num_residues,
                                            connected_domain_ext_msf, connected_domain_int_msf, domains[d].num_residues)
            print(f"Connected domains ({fixed_domain_id} - {d}) ratio = {ratio}")
            if ratio < self.min_ratio:
                print("Ratio below minimum criteria. Break.")
                return False
        print("All ratios met the minimum criteria.")
        return True

    def calc_domain_int_msf(self, domain_residues: int, r: gemmi.SupResult):
        """
        Calculates the domain's internal Mean Square Fluctuation
        :param domain_residues: The specified domain's number of residues
        :param r: The gemmi.SupResult associated with the domain
        :return:
        """
        return (r.rmsd ** 2) * domain_residues

    def calc_domain_ext_msf(self, domain, r, protein_id=1):
        """
        Calculates the domain's external Mean Square Fluctuation. The function first transforms the domain chain using
        the given r, then calculates the displacement vectors between atoms of the original domain chain and the
        transformed domain chain.
        :param domain: A domain of the chain
        :param r: The gemmi.SupResult associated with the domain
        :param protein_id: The ID of the Protein used
        :return:
        """
        protein: Protein = self.protein_1 if protein_id == 1 else self.protein_2
        slide_residues = protein.get_slide_window_result()
        domain_chain: gemmi.Chain = gemmi.Chain(protein.chain_param)
        transformed_chain: gemmi.Chain = gemmi.Chain(protein.chain_param)
        for s in domain.segments:
            for i in range(s[0], s[1]+1):
                domain_chain.add_residue(slide_residues[i])
                transformed_chain.add_residue(slide_residues[i])
        transformed_polymer: gemmi.ResidueSpan = transformed_chain.get_polymer()
        transformed_polymer.transform_pos_and_adp(r.transform)
        ext_msf = 0
        for r in range(len(domain_chain)):
            for a in self.backbone_atoms:
                atom_coords = np.asarray(domain_chain[r].sole_atom(a).pos.tolist())
                transformed_atom_coords = np.asarray(transformed_chain[r].sole_atom(a).pos.tolist())
                disp_vec = (atom_coords - transformed_atom_coords) ** 2
                sum_disp = np.sum(disp_vec)
                ext_msf += sum_disp

        return ext_msf

    def calc_ext_int_ratio(self, fixed_ext, fixed_int, fixed_num_residues,
                           connected_ext, connected_int, connected_num_residues):
        sum_exts = (fixed_ext/fixed_num_residues) + (connected_ext/connected_num_residues)
        sum_ints = (fixed_int/fixed_num_residues) + (connected_int/connected_num_residues)
        ratio = math.sqrt(sum_exts/sum_ints)
        return ratio

    def check_domain_sizes(self, domains):
        return False if any(d.num_residues < self.min_domain_size for d in domains) else True

    def check_cluster_sizes(self, segments: dict):
        for segments_list in segments.values():
            res_count = sum([segment[1] + 1 - segment[0] for segment in segments_list])
            if res_count < self.min_domain_size:
                return False
        return True

    def print_labels(self):
        print("Printing Labels...")
        print(f"Length = {len(self.k_means_results.labels_)}")
        print(f"Labels = {[self.k_means_results.labels_[i] for i in range(len(self.k_means_results.labels_))]}")

    def print_segments(self, segments):
        print("SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS")
        print("Printing segments...")
        seg_count = 0
        res_total = 0
        total_small_segments = 0
        for k, v in segments.items():
            print(f"Cluster {k} - {len(v)} Segments")
            print(f"Values {v}")
            res_count = 0
            seg_count += v.shape[0]
            for i in v:
                # print(self.k_means_results.labels_[i[0]:i[1]+1])
                res_count += i[1] + 1 - i[0]
                res_total += i[1] + 1 - i[0]
            print(f"Segment Res = {res_count}")
            if res_count < self.min_domain_size:
                total_small_segments += 1
        print("---------------------------------------")
        print(f"Seg total = {seg_count}")
        print(f"Res total = {res_total}")
        print(f"Total small segments = {total_small_segments}")
        print("SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS")
        return

    def print_domains(self, domains):
        print("DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD")
        print("Printing domains")
        for d in domains:
            print(d)
        print("DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD")

    def print_coords(self, transformed_1, transformed_2):
        for i in range(len(transformed_1)):
            for a in self.backbone_atoms:
                print(f"==============================================")
                print(f"P1 : {self.protein_1.slide_window_residues[i].sole_atom(a).pos} -> {transformed_1[i].sole_atom(a).pos}")
                print(f"P2 : {self.protein_2.slide_window_residues[i].sole_atom(a).pos} -> {transformed_2[i].sole_atom(a).pos}")
                print(f"==============================================")

    def print(self):
        print("Printing Clustered Result")
        print("=========================================")
        print("Printing Domains")
        for d in self.domains:
            print(d)
        print(f"Fixed domain ID = {self.fixed_domain}")
        print("=========================================")


