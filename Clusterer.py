import math
import numpy as np
import gemmi
import Protein
import DomainBuilder as dom_build
from Domain import Domain
from hkmeans import HKMeans


class Clusterer:
    def __init__(self, params: dict, rotation_vectors: np.array, protein_1: Protein, protein_2: Protein, main_atoms):
        # We want a large enough k so that it doesn't get used. If it does manage to get used, good luck.
        self.max_k = 40
        self.current_k = 2
        self.min_domain_size: int = int(params["domain"])
        self.min_ratio: float = float(params["ratio"])
        self.backbone_atoms = main_atoms
        self.num_atoms = len(self.backbone_atoms)
        self.rotation_vectors: np.array = rotation_vectors
        self.protein_1: Protein = protein_1
        self.protein_2: Protein = protein_2
        self.k_means_results = None
        self.valid_cluster = 0
        self.segments = {}
        self.domains = []
        # The domain ID number determining which domain will be the fixed point
        self.fixed_domain = None
        self.clusterer_status = 1

    def cluster(self):
        fail_limit = 5
        # This condition determines whether a k-value cluster is valid based on number of residues
        valid_cluster_found = False
        # This condition determines whether a valid cluster contains domains that have valid domain pair ratio
        valid_domains_found = False
        while self.current_k < self.max_k:
            print("============================================")
            print(f"\ncurrent_k = {self.current_k}")
            # KMeans the rotation vectors to obtain k number of clusters
            self.k_means_results = self.calc_k_means()
            # self.print_labels()
            # Obtain the segments from the Hartigan KMeans results.
            # The segments' indices are for the slide windowed residues.
            temp_segments, cluster_residues_small = self.determine_segments()
            # self.print_segments(temp_segments)
            # If there is a cluster where its total residue is smaller than min domain size:
            # If there is a previous iteration where valid clusters and valid domain pair ratios are found, clustering
            # can be halted
            if cluster_residues_small and valid_domains_found:
                print("Found cluster with number of residues smaller than minimum domain size. Clustering halted. ")
                print(f"Final clusterer to use: {self.valid_cluster}")
                self.clusterer_status = 0
                return
            # If there are no previous iterations where domains have been found yet, skip to the next k value
            if cluster_residues_small and not valid_cluster_found:
                print("Found cluster with number of residues smaller than minimum domain size. Going to next k value")
                self.current_k += 1
                if self.current_k > fail_limit:
                    print("Too many fails. Increasing window size by 2")
                    self.clusterer_status = -1
                    return
                continue
            # If there was a previous iteration where a cluster had domains but the domain pair ratios are invalid,
            # reset the clustering but this time increase the window size by 2.
            if cluster_residues_small and valid_cluster_found:
                print("Found cluster with number of residues smaller than minimum domain size. Valid cluster previously "
                      "found but no domains. Increasing window size by 2")
                self.current_k += 1
                self.clusterer_status = -1
                return
            valid_cluster_found = True
            self.print_segments(temp_segments)
            temp_domains, cluster_break = dom_build.build_domains(self.protein_1.slide_window_residues,
                                                                  self.protein_2.slide_window_residues,
                                                                  temp_segments, self.min_domain_size)
            # The first time that a k value produces valid domains for all clusters, set it to true
            if cluster_break and valid_domains_found:
                print("All domains found are smaller than minimum domain size. Clustering halted.")
                self.clusterer_status = 0
                return
            elif cluster_break and not valid_cluster_found:
                print("All domains in the cluster are less than minimum domain size")
                self.current_k += 1
                continue
            # self.print_domains(temp_domains_1, current_k)
            temp_fixed_domain_id = self.check_domain_connectivity(temp_domains)
            # Perform mass-weighted whole-protein best fit to obtain a new set of Protein 2 coordinates after
            # superimposing Protein 2 slide-window chain onto Protein 1 slide-window chain
            # _, r2, _, transformed_2_on_1 = self.mass_weighted_whole_protein_fit(temp_domains_1)

            ratio_not_met = False
            for domain in temp_domains:
                if domain.domain_id == temp_fixed_domain_id:
                    continue
                transformed_2_domains_on_1 = self.mass_weighted_fit(domain, temp_domains[temp_fixed_domain_id])
                ratio_met = self.check_ratios(transformed_2_domains_on_1, domain, temp_domains[temp_fixed_domain_id])
                if not ratio_met:
                    ratio_not_met = True
                    print("Ratio not met")
                    break
            if ratio_not_met:
                self.current_k += 1
                continue
            else:
                self.valid_cluster = self.current_k
                self.domains = temp_domains
                self.fixed_domain = temp_fixed_domain_id
                self.segments = temp_segments
                valid_domains_found = True
            self.current_k += 1

    def calc_k_means(self):
        """
        Performs Hartigan-Wong KMeans on the rotation vectors
        :return: KMeans results
        """
        k_means = HKMeans(n_clusters=self.current_k, n_init=15, max_iter=300).fit_predict(self.rotation_vectors)
        return k_means

    def determine_segments(self):
        """
        Finds segments in the array of K-Means cluster IDs where the segments contain the same cluster IDs in a row.
        :return:    A dictionary where the keys are the cluster IDs and the values are a list of tuples (start, end)
                    where start is the index in k_mean_labels where the segment starts and end is the index in
                    k_mean_labels where the segment ends.
        """
        # print(len(k_mean_labels))
        # Set the first label as the label checker
        current_element_to_check = self.k_means_results[0]
        start_index = 0
        end_index = 0

        # This line of code to declare and initialise a dictionary of lists does not work as it will end up
        # appending to all lists.
        # segment_indices = dict.fromkeys(range(0, k), [])

        # This is the correct way to initialise a dictionary of lists
        segment_indices = {key: np.array([[]], dtype=int) for key in range(self.current_k)}

        # Iterate through each label
        for i in range(len(self.k_means_results)):
            # When the label is not equal to the checker. It means the segment ends there and the segment's start and
            # end indices are obtained and stored.
            if self.k_means_results[i] != current_element_to_check:
                temp = np.append(segment_indices[current_element_to_check], [[start_index, end_index]],
                                 axis=0 if segment_indices[current_element_to_check].shape[1] > 0 else 1)
                segment_indices[current_element_to_check] = temp
                current_element_to_check = self.k_means_results[i]
                start_index = i
            end_index = i
            if i == len(self.k_means_results) - 1:
                temp = np.append(segment_indices[current_element_to_check], [[start_index, end_index]],
                                 axis=0 if segment_indices[current_element_to_check].shape[1] > 0 else 1)
                segment_indices[current_element_to_check] = temp
                # segment_indices[current_element_to_check].append((start_index, i))

        cluster_residues_too_little = False
        # Checks whether each cluster's total residues have minimum domain size. If a cluster is less than minimum
        # domain size, it means that the clustering will stop.
        # print(f"Values {segment_indices.values()}")
        for indices in segment_indices.values():
            # print(indices)
            num_residues = indices[:, 1] + 1 - indices[:, 0]
            # print(f"Num res = {num_residues}")
            # print(f"Sum = {sum(num_residues)}")
            if sum(num_residues) < 20:
                cluster_residues_too_little = True
                break

        return segment_indices, cluster_residues_too_little

    def check_domain_connectivity(self, domains):
        """
        Takes the domains of one of the protein conformation and determines which domain has the most number of domains
        connected to it.
        :return chosen_domain: The ID of the domain with the most connected domains
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

    def mass_weighted_fit(self, connected_domain: Domain, fixed_domain: Domain):
        """
        Performs a mass-weighted protein best fit on a domain pair to get a new set of coordinates of a domain pair.
        :param connected_domain: The domain connected to the fixed domain
        :param fixed_domain: The fixed domain
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
        coords_1 = []
        coords_2 = []
        weights = []

        connected_domain_num_residues = connected_domain.num_residues * len(self.backbone_atoms)
        fixed_domain_num_residues = fixed_domain.num_residues * len(self.backbone_atoms)

        for s in connected_domain.segments:
            for i in range(s[0], s[1] + 1):
                for a in self.backbone_atoms:
                    coords_1.append(slide_window_1[i].sole_atom(a).pos)
                    coords_2.append(slide_window_2[i].sole_atom(a).pos)
                    weights.append(1 / connected_domain_num_residues)

        for s in fixed_domain.segments:
            for i in range(s[0], s[1] + 1):
                for a in self.backbone_atoms:
                    coords_1.append(slide_window_1[i].sole_atom(a).pos)
                    coords_2.append(slide_window_2[i].sole_atom(a).pos)
                    weights.append(1 / fixed_domain_num_residues)

        # The superposition result of Protein 1 fitting onto Protein 2
        # r_1: gemmi.SupResult = gemmi.superpose_positions(coords_2, coords_1, weights)
        # The superposition result of Protein 2 fitting onto Protein 1
        r_2: gemmi.SupResult = gemmi.superpose_positions(coords_1, coords_2, weights)
        # slide_window_1.transform_pos_and_adp(r_1.transform)
        slide_window_2.transform_pos_and_adp(r_2.transform)
        # print(f"Count = {r_2.count}")
        # print(f"Center 1 = {r_2.center1}")
        # print(f"Center 2 = {r_2.center2}")
        # print(f"RMSD = {r_2.rmsd}")
        # print(f"Rotation = {r_2.transform.mat}")
        # print(f"Translation = {r_2.transform.vec}")
        # return r_1, r_2, slide_window_1, slide_window_2
        return slide_window_2

    def check_ratios(self, transformed_protein, connected_domain: Domain, fixed_domain: Domain):
        """
        Checks the ratio of internal and external domain movements of each fixed-connected domain pairs.
        :param transformed_protein: The mass-weight fitted Protein slide window chain onto the other Protein.
        :param connected_domain: The domain connected to the fixed domain
        :param fixed_domain: The fixed domain
        :return:
        """
        fixed_domain_r: gemmi.SupResult = self.superimpose_domain(transformed_protein, fixed_domain)
        fixed_domain_num_atoms = fixed_domain.num_residues * len(self.backbone_atoms)
        # fixed_domain_num_atoms = domains[fixed_domain_id].num_residues
        fixed_domain_int_msf = self.calc_domain_int_msf(fixed_domain_num_atoms, fixed_domain_r)
        fixed_domain_ext_msf = self.calc_domain_ext_msf(fixed_domain, fixed_domain_r)

        connected_domain_r: gemmi.SupResult = self.superimpose_domain(transformed_protein, connected_domain)
        domain_num_atoms = connected_domain.num_residues * len(self.backbone_atoms)
        # domain_num_atoms = domains[d].num_residues
        # connected_domain_int_msf = self.calc_domain_int_msf(domains[d].num_residues, results[d])
        connected_domain_int_msf = self.calc_domain_int_msf(domain_num_atoms, connected_domain_r)
        connected_domain_ext_msf = self.calc_domain_ext_msf(connected_domain, connected_domain_r)
        ratio = self.calc_ext_int_ratio(fixed_domain_ext_msf, fixed_domain_int_msf, fixed_domain_num_atoms,
                                        connected_domain_ext_msf, connected_domain_int_msf, domain_num_atoms)
        print(f"Connected domains ({fixed_domain.domain_id} - {connected_domain.domain_id}) ratio = {ratio}")
        if ratio < self.min_ratio:
            print("Ratio below minimum criteria. Break.")
            return False
        else:
            return True

    def superimpose_domain(self, residue_span, domain: Domain, protein=1):
        """
        Superimposes each Protein domain's original residue atoms onto the transformed residue atoms
        :param domain: Domain object
        :param residue_span: A ResidueSpan object
        :param protein: The ID of the Protein to be used to superimpose onto residue_span
        :return results: A list of gemmi.SupResult objects
        """
        # Get the original chain from the Protein that will be used to superimpose onto the target
        protein: Protein = self.protein_1 if protein == 1 else self.protein_2
        fitting_chain: gemmi.ResidueSpan = protein.get_slide_window_result()
        # Original residue atoms as the fitting
        fitting_coords = []
        # The transformed residue atoms as the target
        target_coords = []
        for s in domain.segments:
            for si in range(s[0], s[1] + 1):
                for a in self.backbone_atoms:
                    fitting_coords.append(fitting_chain[si].sole_atom(a).pos)
                    target_coords.append(residue_span[si].sole_atom(a).pos)
        # Superpose the domain residue atoms onto the target residue atoms
        r = gemmi.superpose_positions(target_coords, fitting_coords)
        return r

    def calc_domain_int_msf(self, domain_atoms: int, r: gemmi.SupResult):
        """
        Calculates the domain's internal Mean Square Fluctuation
        :param domain_atoms: The specified domain's number of atoms
        :param r: The gemmi.SupResult associated with the domain
        :return:
        """
        return (r.rmsd ** 2) * domain_atoms

    def calc_domain_ext_msf(self, domain, r, protein_id=1):
        """
        Calculates the domain's external Mean Square Fluctuation. The function first transforms the domain chain using
        the given r, then calculates the displacement vectors between atoms of the original domain chain and the
        transformed domain chain.
        :param domain: A Domain object
        :param r: The gemmi.SupResult associated with the domain
        :param protein_id: The ID of the Protein used
        :return:
        """
        protein: Protein = self.protein_1 if protein_id == 1 else self.protein_2
        slide_residues = protein.get_slide_window_result()
        original_chain: gemmi.Chain = gemmi.Chain(protein.chain_param)
        transformed_chain: gemmi.Chain = gemmi.Chain(protein.chain_param)
        for s in domain.segments:
            for i in range(s[0], s[1]+1):
                original_chain.add_residue(slide_residues[i])
                transformed_chain.add_residue(slide_residues[i])
        transformed_polymer: gemmi.ResidueSpan = transformed_chain.get_polymer()
        transformed_polymer.transform_pos_and_adp(r.transform)
        ext_msf = 0
        for r in range(len(original_chain)):
            for a in self.backbone_atoms:
                atom_coords = np.asarray(original_chain[r].sole_atom(a).pos.tolist())
                transformed_atom_coords = np.asarray(transformed_chain[r].sole_atom(a).pos.tolist())
                disp_vec = atom_coords - transformed_atom_coords
                sum_disp = np.sum(disp_vec ** 2)
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

    def renumber_domains(self):
        pass

    def print_labels(self):
        print("Printing Labels...")
        print(f"Length = {len(self.k_means_results)}")
        curr_element = self.k_means_results[0]
        list_of_same_labels = []
        for i in range(len(self.k_means_results)):
            if self.k_means_results[i] != curr_element:
                print(list_of_same_labels)
                list_of_same_labels = [self.k_means_results[i]]
                curr_element = self.k_means_results[i]
                continue
            list_of_same_labels.append(self.k_means_results[i])
        print(list_of_same_labels)

    def print_segments(self, segments):
        print("SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS")
        print("Printing segments...")
        seg_count = 0
        res_total = 0
        total_small_segments = 0
        for k, v in segments.items():
            print(f"Cluster {k} : {len(v)} Segments")
            print(f"Values {v}")
            res_count = 0
            seg_count += v.shape[0]
            for i in v:
                res_count += i[1] + 1 - i[0]
                res_total += i[1] + 1 - i[0]
            print(f"Segment Res = {res_count}")
            if res_count < self.min_domain_size:
                total_small_segments += 1
        print("---------------------------------------")
        print(f"Seg total = {seg_count}")
        print(f"Res total = {res_total}")
        print("SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS")
        return

    def print_coords(self, transformed_1, transformed_2):
        for i in range(len(transformed_1)):
            for a in self.backbone_atoms:
                print(f"==============================================")
                print(f"P1 : {self.protein_1.slide_window_residues[i].sole_atom(a).pos} -> {transformed_1[i].sole_atom(a).pos}")
                print(f"P2 : {self.protein_2.slide_window_residues[i].sole_atom(a).pos} -> {transformed_2[i].sole_atom(a).pos}")
                print(f"==============================================")

