import math
import numpy as np
import gemmi
import difflib
import FileMngr
from Clusterer import Clusterer
from Protein import Protein
from scipy.spatial.transform import Rotation


class Engine:
    def __init__(self, first_pdb: str, second_pdb: str, params=None):
        # Default parameters to be used if not given input
        if params is None:
            self.parameters = {"chain1id": "A",
                               "chain2id": "A",
                               "window": "5",
                               "domain": "20",
                               "ratio": "1.0",
                               "atoms": "backbone"}
        else:
            self.parameters = params
        # Initialise proteins
        self.protein_1: Protein = Protein(first_pdb, self.parameters["chain1id"], self.parameters["atoms"])
        self.protein_2: Protein = Protein(second_pdb, self.parameters["chain2id"], self.parameters["atoms"])
        self.fitting_protein_1 = None
        self.main_atoms = []
        if self.parameters["atoms"] == "backbone":
            self.main_atoms = ["N", "CA", "C"]
        elif self.parameters["atoms"] == "ca":
            self.main_atoms = ["CA"]
        # A object containing superimposition information between the proteins
        self.chain_superimpose_result = None
        # List of arrays of gemmi.SupResult objects containing superimposition information between each residue
        self.slide_window_superimpose_results = []
        # Array of translation vectors
        self.translation_vecs = None
        # Array of (3 x 3) rotation matrices
        self.rotation_mats = None
        # List of rotation vectors created from rotation_mats
        self.rotation_vecs = None
        # List of unit vectors of the rotation vectors
        self.unit_vectors = None
        # List of angles
        self.angles = None
        # Sum of all angles
        self.angles_sum = 0.0
        # Scaling of rotation vectors
        self.scaling = 0.0
        # Clusterer object
        self.clusterer = None

    def run(self):
        """
        Runs the entire program. To make it simple, all operations happen based on Protein 2 conforming to Protein 1.
        :return:
        """
        if self.check_chain_compatibility_and_length():
            # Superimposes Protein 2 onto Protein 1 so they are in the same "space". The superimposed protein will be
            # called Protein S.
            self.superimpose_chains()
            # self.print_chains_superimposed_result()
            # Slides a window over Protein 1's residues and Protein S's residues and superimposes them in each window.
            self.sliding_window_superimpose_residues()
            # self.sliding_window_superimpose_residues(fitted_protein_polymer)
            # Obtain the rotation matrices from the superposition results
            self.get_transformations()
            # Convert rotation matrices to rotation vectors
            self.convert_rot_mats()
            # self.print_rotation_vectors()
            FileMngr.write_rotation_vec_to_pdb(self.protein_1.id, self.protein_1.slide_window_residues,
                                               self.protein_1.slide_window_residues_indices, self.rotation_vecs)
            # self.print_unit_vectors()
            cluster_status = 1
            while cluster_status != 0:
                self.clusterer: Clusterer = Clusterer(self.parameters, self.rotation_vecs,
                                                      self.protein_1, self.protein_2, self.main_atoms)
                cluster_status = self.clusterer.cluster()
                if cluster_status == -1:
                    self.parameters["window"] = self.parameters["window"] + 2
            self.determine_domains_screw_axis()

        else:
            print("Unable to compare sequences.")
            return False
        return True

    def check_chain_compatibility_and_length(self):
        """
        Checks whether the residues sequence of Protein 1 and Protein 2 are similar. The function also checks whether
        the chains are of the same length. If they are not but the difference is small, remove some residues.
        Note: Need to find a way to check which residues to remove.
        :return: Boolean on whether the chains are compatible
        """
        # residues_1 = [res.name for res in self.protein_1.residue_span]
        # residues_2 = [res.name for res in self.protein_2.residue_span]
        # sm = difflib.SequenceMatcher(None, residues_1, residues_2)
        # if self.protein_1.chain_atoms.shape[0] == self.protein_2.chain_atoms.shape[0] and sm.ratio() > 0.4:

        # If the protein chains have identical lengths, then nothing else needs to be done.
        if self.protein_1.chain_atoms.shape[0] == self.protein_2.chain_atoms.shape[0]:
            window_size = int(self.parameters["window"])
            window_mid_index = start_index = (window_size - 1) // 2
            protein_1_residues_length = len(self.protein_1.get_polymer())
            # The index of the middle residue in the last sliding window
            final_index = protein_1_residues_length + 1 - (window_size - window_mid_index)
            self.protein_1.slide_window_residues_indices = (start_index, final_index)
            self.protein_2.slide_window_residues_indices = (start_index, final_index)
            return True

        if abs(self.protein_2.chain_atoms.shape[0] - self.protein_1.chain_atoms.shape[0]) > 5:
            print("Chains are too different. AKA: One of them is too long.")
            return False
        else:
            window_size = int(self.parameters["window"])
            window_mid_index = start_index = (window_size - 1) // 2
            if self.protein_1.chain_atoms.shape[0] < self.protein_2.chain_atoms.shape[0]:
                protein_1_residues_length = len(self.protein_1.get_polymer())
                # The index of the middle residue in the last sliding window
                final_index = protein_1_residues_length + 1 - (window_size - window_mid_index)
            else:
                protein_2_residues_length = len(self.protein_2.get_polymer())
                # The index of the middle residue in the last sliding window
                final_index = protein_2_residues_length + 1 - (window_size - window_mid_index)
            self.protein_1.slide_window_residues_indices = (start_index, final_index)
            self.protein_2.slide_window_residues_indices = (start_index, final_index)
            return True

    def superimpose_chains(self):
        """
        Superimposes the entire chain of Protein 2 unto Protein 1 using the backbone atoms.
        :return: Protein 2 chain after transformation
        """
        ptype = self.protein_1.residue_span.check_polymer_type()

        protein_1 = self.protein_1.get_polymer()
        coords_1 = [r.sole_atom(a).pos for r in protein_1 for a in self.main_atoms]
        protein_2 = self.protein_2.get_polymer()
        coords_2 = [r.sole_atom(a).pos for r in protein_2 for a in self.main_atoms]

        # fit_1_to_2_result: gemmi.SupResult = gemmi.calculate_superposition(self.protein_2.get_slide_window_result(),
        #                                                                    self.protein_1.get_slide_window_result(),
        #                                                                    ptype, gemmi.SupSelect.MainChain)
        fit_1_to_2_result: gemmi.SupResult = gemmi.superpose_positions(coords_2, coords_1)
        self.chain_superimpose_result = fit_1_to_2_result
        self.fitting_protein_1: gemmi.ResidueSpan = self.protein_1.get_polymer()
        self.fitting_protein_1.transform_pos_and_adp(fit_1_to_2_result.transform)
        # print(f"Count = {fit_1_to_2_result.count}")
        # print(f"Center 1 = {fit_1_to_2_result.center1}")
        # print(f"Center 2 = {fit_1_to_2_result.center2}")
        # print(f"RMSD = {fit_1_to_2_result.rmsd}")
        # print(f"Rotation = {fit_1_to_2_result.transform.mat}")
        # print(f"Translation = {fit_1_to_2_result.transform.vec}")

    def sliding_window_superimpose_residues(self):
        """
        Slides an array window of a specified size over both protein chains to superimpose the backbone atoms of the
        residues inside the window and getting indices of each residue located in the middle of each window.
        :return:
        """
        residues_length = len(self.fitting_protein_1)  # The number of residues in Protein S
        window_size = int(self.parameters["window"])
        target_protein_polymer = self.protein_2.get_polymer()  # The protein that Protein S will superimpose onto.

        # Only used if the window size ends up becoming larger or equal to the length of the chain. Hopefully, this
        # condition isn't met since the length of a protein chain is quite long.
        if residues_length <= window_size:
            target_protein_polymer_atoms_pos = []
            fitting_protein_polymer_atoms_pos = []
            # The middle index of the large sliding window
            residues_mid_index = (residues_length - 1) // 2
            indices = (residues_mid_index, residues_mid_index + 1)
            self.protein_1.slide_window_residues_indices = indices
            self.protein_2.slide_window_residues_indices = indices
            self.get_utilised_residues(indices)
            for r in range(len(target_protein_polymer)):
                for a in self.main_atoms:
                    target_protein_polymer_atoms_pos.append(target_protein_polymer[r].sole_atom(a).pos)
                    fitting_protein_polymer_atoms_pos.append(self.fitting_protein_1[r].sole_atom(a).pos)
            self.slide_window_superimpose_results.append(gemmi.superpose_positions(target_protein_polymer_atoms_pos,
                                                                                   fitting_protein_polymer_atoms_pos))
        else:
            self.get_utilised_residues(self.protein_1.slide_window_residues_indices)
            # For each residue in the 2 proteins, get the backbone atoms
            # Instead of going from the first residue to the last residue of the last sliding window, use the already
            # calculated slide_window_residues_indices.
            if window_size % 2 == 1:
                first_half = second_half = window_size // 2
            else:
                first_half = (window_size / 2) - 1
                second_half = window_size / 2
            for r in range(self.protein_1.slide_window_residues_indices[0], self.protein_1.slide_window_residues_indices[1]):
                target_protein_polymer_atoms_pos = []
                fitting_protein_polymer_atoms_pos = []
                residue_start_index = r - first_half
                residue_end_index = r + second_half + 1
                for i in range(residue_start_index, residue_end_index):
                    for a in self.main_atoms:
                        target_protein_polymer_atoms_pos.append(target_protein_polymer[i].sole_atom(a).pos)
                        fitting_protein_polymer_atoms_pos.append(self.fitting_protein_1[i].sole_atom(a).pos)
                # Superimpose and append the result to list
                self.slide_window_superimpose_results.append(gemmi.superpose_positions(target_protein_polymer_atoms_pos,
                                                                                       fitting_protein_polymer_atoms_pos))

    def get_utilised_residues(self, indices):
        """
        Gets all residues that were used in the middle of each sliding window
        :param indices:
        :return:
        """
        chain_1 = self.protein_1.get_polymer()
        chain_2 = self.protein_2.get_polymer()
        for i in range(indices[0], indices[1]):
            self.protein_1.slide_window_residues.append(chain_1[i])
            self.protein_2.slide_window_residues.append(chain_2[i])

    def get_transformations(self):
        """
        residues_superimpose_results is a list of gemmi.SupResult objects containing information of the superimposition
        between each sliding window of the Proteins. This function is used to extract the numerical data of the objects
        for KMeans clustering. Specifically, the rotation matrix.
        :return:
        """
        # Initialise the numpy arrays
        # Array of rotation matrices
        self.rotation_mats = np.empty(shape=[len(self.slide_window_superimpose_results), 3, 3])
        # Array of translation vectors
        self.translation_vecs = np.empty(shape=[len(self.slide_window_superimpose_results), 3])
        for i in range(len(self.slide_window_superimpose_results)):
            rot_mat: gemmi.Mat33 = self.slide_window_superimpose_results[i].transform.mat
            trans_vec: gemmi.Vec3 = self.slide_window_superimpose_results[i].transform.vec
            self.rotation_mats[i] = np.asarray(rot_mat.tolist())
            self.translation_vecs[i] = np.asarray(trans_vec.tolist())

    def convert_rot_mats(self):
        """
        Convert rotation matrices to rotation vectors, angles in degrees, and unit rotation vectors
        :return:
        """
        self.rotation_vecs = Rotation.from_matrix(self.rotation_mats).as_rotvec(degrees=True)
        self.angles = np.array([np.linalg.norm(i) for i in self.rotation_vecs])
        self.unit_vectors = self.rotation_vecs / np.linalg.norm(self.rotation_vecs)
        self.angles_sum = np.sum(self.angles)
        self.scaling = 20 * (self.angles_sum / (self.protein_1.slide_window_residues_indices[1] - self.protein_1.slide_window_residues_indices[0]))

    def determine_domains_screw_axis(self):
        fixed_domain_transformations = self.clusterer.rs
        original_slide_chain_2: gemmi.ResidueSpan = self.protein_2.get_slide_window_result()
        transformed_slide_chain_1: gemmi.ResidueSpan = self.protein_1.get_slide_window_result()
        transformed_slide_chain_1.transform_pos_and_adp(fixed_domain_transformations.transform)

        domain_disp_vecs = []

        for domain in self.clusterer.domains:
            original_domain_chain: gemmi.Chain = gemmi.Chain(self.protein_2.chain_param)
            transformed_domain_chain: gemmi.Chain = gemmi.Chain(self.protein_1.chain_param)
            if domain.domain_id == self.clusterer.fixed_domain:
                continue
            original_atoms = []
            transformed_atoms = []
            for segment in domain.segments:
                for s in segment:
                    for i in range(s[0], s[1] + 1):
                        original_domain_chain.add_residue(original_slide_chain_2[i])
                        transformed_domain_chain.add_residue(original_slide_chain_2[i])
                        for a in self.main_atoms:
                            original_atoms.append(original_slide_chain_2[i].sole_atom(a).pos)
                            transformed_atoms.append(transformed_slide_chain_1[i].sole_atom(a).pos)

            r: gemmi.SupResult = gemmi.superpose_positions(original_atoms, transformed_atoms)
            original_domain_residues = original_domain_chain.get_polymer()
            transformed_domain_residues: gemmi.ResidueSpan = transformed_domain_chain.get_polymer()
            transformed_domain_residues.transform_pos_and_adp(r.transform)
            disp_vecs = np.array((len(original_domain_residues), 3))
            for r in range(len(original_domain_residues)):
                for a in range(len(self.main_atoms)):
                    atom_coords = np.asarray(original_domain_residues[r].sole_atom(self.main_atoms[a]).pos.tolist())
                    transformed_atom_coords = np.asarray(transformed_domain_residues[r].sole_atom(self.main_atoms[a]).pos.tolist())
                    disp_vecs[(r*3)+a] = atom_coords - transformed_atom_coords

            norm_disp_vecs = disp_vecs / np.linalg.norm(disp_vecs)
            domain_disp_vecs.append(norm_disp_vecs)




        return

    def determine_external_component_motion(self):

        return

    def print_chains_superimposed_result(self):
        # A gemmi.SupResult object containing superimposition information between 2 chains
        print(f"RMSD =                  {self.chain_superimpose_result.rmsd}")
        print(f"Count =                 {self.chain_superimpose_result.count}")
        print(f"Center 1 =              {self.chain_superimpose_result.center1}")
        print(f"Center 2 =              {self.chain_superimpose_result.center2}")
        print(f"Translation Vector =    {self.chain_superimpose_result.transform.vec}")
        print(f"Rotation Matrix =       {self.chain_superimpose_result.transform.mat}")

    def print_slide_window_superimpose_results(self, n=None):
        if n is None or n > len(self.slide_window_superimpose_results):
            n = len(self.slide_window_superimpose_results)
        print(f"slide_window_superimpose_result size = {len(self.slide_window_superimpose_results)}")
        for i in range(n):
            item: gemmi.SupResult = self.slide_window_superimpose_results[i]
            print(f"RMSD =                  {item.rmsd}")
            print(f"Count =                 {item.count}")
            print(f"Center 1 =              {item.center1}")
            print(f"Center 2 =              {item.center2}")
            print(f"Translation Vector =    {item.transform.vec}")
            print(f"Rotation Matrix =       {item.transform.mat}")

    def print_slide_window_residue_indices(self):
        print(f"slide_window_residue_indices = {self.protein_1.slide_window_residues_indices}")

    def print_rotation_matrices(self, n=None):
        if n is None or n > self.rotation_mats.shape[0]:
            n = self.rotation_mats.shape[0]
        print(f"rotation_mats shape = {self.rotation_mats.shape}")
        print(f"rotation_mats[0:{n}] = {self.rotation_mats[0:n]}")

    def print_rotation_vectors(self, n=None):
        print(f"rotation_vecs shape = {self.rotation_vecs.shape}")
        if n is None or n > self.rotation_vecs.shape[0]:
            n = self.rotation_vecs.shape[0]
        for i in range(n):
            print(f"{self.protein_1.slide_window_residues[i]} - [{self.rotation_vecs[i][0]}, {self.rotation_vecs[i][1]}, {self.rotation_vecs[i][2]}]")

    def print_unit_vectors(self, n=None):
        print(f"rotation_vecs shape = {self.unit_vectors.shape}")
        if n is None or n > self.unit_vectors.shape[0]:
            n = self.unit_vectors.shape[0]
        for i in range(n):
            print(
                f"{self.protein_1.slide_window_residues[i]} - [{self.unit_vectors[i][0]}, {self.unit_vectors[i][1]}, {self.unit_vectors[i][2]}]")

    # def print_angles(self, n=None):
    #     if n is None or n > self.angles.shape[0]:
    #         n = self.angles.shape[0]
    #     print(f"angles shape = {self.angles.shape}")
    #     print(f"angles[0:{n} = {self.angles[0:n]}")

    #
    # def print_chains_superimposed_result(self):
    #     print(math.sqrt(sum([r.rmsd for r in self.chain_superimpose_result])))
    #     # for r in self.chain_superimpose_result:
    #     #     print(f"RMSD =                  {r.rmsd}")
    #     #     print(f"Count =                 {r.count}")
    #     #     print(f"Center 1 =              {r.center1}")
    #     #     print(f"Center 2 =              {r.center2}")
    #     #     print(f"Translation Vector =    {r.transform.vec}")
    #     #     print(f"Rotation Matrix =       {r.transform.mat}")

