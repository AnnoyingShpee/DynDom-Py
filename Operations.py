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
        self.fitting_protein_2 = None
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
            self.clusterer: Clusterer = Clusterer(self.parameters, self.rotation_vecs,
                                                  self.protein_1, self.protein_2, self.main_atoms)
            self.clusterer.cluster()
            # self.clusterer.print()

        else:
            print("Unable to compare sequences.")
            return False
        return True

    def check_chain_compatibility_and_length(self):
        """
        Checks whether the residues sequence of Protein 1 and Protein 2 are similar
        :return: The percentage of similarity of the 2 protein residue sequences (Between 0 and 1)
        """
        residues_1 = [res.name for res in self.protein_1.residue_span]
        residues_2 = [res.name for res in self.protein_2.residue_span]
        sm = difflib.SequenceMatcher(None, residues_1, residues_2)
        if self.protein_1.chain_atoms.shape[0] == self.protein_2.chain_atoms.shape[0] and sm.ratio() > 0.4:
            return True
        else:
            return False

    def superimpose_chains(self):
        """
        Superimposes the entire chain of Protein 2 unto Protein 1 using the backbone atoms.
        :return: Protein 2 chain after transformation
        """
        ptype = self.protein_1.residue_span.check_polymer_type()
        fit_2_to_1_result: gemmi.SupResult = gemmi.calculate_superposition(self.protein_1.residue_span,
                                                                           self.protein_2.residue_span,
                                                                           ptype, gemmi.SupSelect.MainChain)
        self.chain_superimpose_result = fit_2_to_1_result
        self.fitting_protein_2: gemmi.ResidueSpan = self.protein_2.get_polymer()
        self.fitting_protein_2.transform_pos_and_adp(fit_2_to_1_result.transform)

    def sliding_window_superimpose_residues(self):
        """
        Slides an array window of a specified size over both protein chains to superimpose the backbone atoms of the
        residues inside the window and getting indices of each residue located in the middle of each window.
        :return:
        """
        residues_length = len(self.fitting_protein_2)  # The number of residues in Protein S
        window_size = int(self.parameters["window"])
        # The index of the middle residue of the window
        window_mid_index = start_index = (window_size - 1) // 2
        target_protein_polymer = self.protein_1.get_polymer()  # The protein that Protein S will superimpose onto.

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
                    fitting_protein_polymer_atoms_pos.append(self.fitting_protein_2[r].sole_atom(a).pos)
            self.slide_window_superimpose_results.append(gemmi.superpose_positions(target_protein_polymer_atoms_pos,
                                                                                   fitting_protein_polymer_atoms_pos))
        else:
            # The index of the middle residue in the last sliding window
            final_index = residues_length + 1 - (window_size - window_mid_index)
            indices = (start_index, final_index)
            self.protein_1.slide_window_residues_indices = indices
            self.protein_2.slide_window_residues_indices = indices
            self.get_utilised_residues(indices)
            # For each residue in the 2 proteins, get the backbone atoms
            # Start from the first residue to the last residue of the last sliding window. This means a few residues
            # at the end of the chain will not be used.
            for r in range(residues_length - window_size + 1):
                target_protein_polymer_atoms_pos = []
                fitting_protein_polymer_atoms_pos = []
                residue_end_index = r + window_size
                for i in range(r, residue_end_index):
                    for a in self.main_atoms:
                        target_protein_polymer_atoms_pos.append(target_protein_polymer[i].sole_atom(a).pos)
                        fitting_protein_polymer_atoms_pos.append(self.fitting_protein_2[i].sole_atom(a).pos)
                # Superimpose and append the result to list
                self.slide_window_superimpose_results.append(gemmi.superpose_positions(target_protein_polymer_atoms_pos,
                                                                                       fitting_protein_polymer_atoms_pos))

    # def sliding_window_superimpose_residues(self, target_protein_polymer: gemmi.ResidueSpan, fitted_protein_id=2):
    #     residues_length = len(target_protein_polymer)
    #     window_size = int(self.parameters["window"])
    #     window_mid_index = start_index = (window_size - 1) // 2
    #     fitting_protein: Protein = self.protein_1 if fitted_protein_id == 2 else self.protein_2
    #     fitting_protein_polymer = fitting_protein.get_polymer()
    #
    #     if residues_length <= window_size:
    #         residues_mid_index = (residues_length - 1) // 2
    #         indices = (residues_mid_index, residues_mid_index + 1)
    #         self.protein_1.slide_window_residues_indices = indices
    #         self.protein_2.slide_window_residues_indices = indices
    #         self.get_utilised_residues(indices)
    #         ptype = target_protein_polymer.check_polymer_type()
    #         self.slide_window_superimpose_results.append(gemmi.calculate_superposition(target_protein_polymer,
    #                                                                                    fitting_protein_polymer,
    #                                                                                    ptype,
    #                                                                                    gemmi.SupSelect.MainChain))
    #
    #     else:
    #         final_index = residues_length + 1 - (window_size - window_mid_index)
    #         indices = (start_index, final_index)
    #         self.protein_1.slide_window_residues_indices = indices
    #         self.protein_2.slide_window_residues_indices = indices
    #         self.get_utilised_residues(indices)
    #         for r in range(residues_length - window_size + 1):
    #             target_chain: gemmi.Chain = gemmi.Chain("A")
    #             fitting_chain: gemmi.Chain = gemmi.Chain("A")
    #             residue_end_index = r + window_size
    #             for i in range(r, residue_end_index):
    #                 target_chain.add_residue(target_protein_polymer[i])
    #                 fitting_chain.add_residue(fitting_protein_polymer[i])
    #             target_chain_polymer: gemmi.ResidueSpan = target_chain.get_polymer()
    #             fitting_chain_polymer: gemmi.ResidueSpan = fitting_chain.get_polymer()
    #             ptype = target_chain_polymer.check_polymer_type()
    #             self.slide_window_superimpose_results.append(gemmi.calculate_superposition(target_chain_polymer,
    #                                                                                        fitting_chain_polymer,
    #                                                                                        ptype,
    #                                                                                        gemmi.SupSelect.MainChain))

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

    def determine_screw_axis(self):
        results = self.clusterer.rs
        r = results[self.clusterer.fixed_domain]
        original_slide_chain_1: gemmi.ResidueSpan = self.protein_1.get_slide_window_result()
        transformed_slide_chain_2: gemmi.ResidueSpan = self.protein_2.get_slide_window_result()
        transformed_slide_chain_2.transform_pos_and_adp(r.transform)
        for domain in self.clusterer.domains:
            if domain.domain_id == self.clusterer.fixed_domain:
                continue
            original_atoms = []
            transformed_atoms = []
            for segment in domain.segments:
                for s in segment:
                    for i in range(s[0], s[1] + 1):
                        for a in self.clusterer.backbone_atoms:
                            original_atoms.append(original_slide_chain_1[i].sole_atom(a).pos)
                            transformed_atoms.append(transformed_slide_chain_2[i].sole_atom(a).pos)


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

