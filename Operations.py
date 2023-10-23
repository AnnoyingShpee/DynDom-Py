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
        # A object containing superimposition information between the proteins
        self.chain_superimpose_result = None
        # List of residues located in the middle of each sliding window
        self.slide_window_residue_indices = ()
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
        # Clusterer object
        self.clusterer = None
        # self.slide_window_result_1, self.slide_window_result_2 = np.array([]), np.array([])

    def run(self):
        if self.check_chain_compatibility() > 0.4 and self.check_chain_lengths():
            # self.superimpose_chains()
            # self.print_chains_superimposed_result()
            self.sliding_window_superimpose_residues()
            # self.print_slide_window_residue_indices()
            # self.print_slide_window_superimpose_results(5)
            # Obtain the rotation matrices from the superposition results
            self.get_transformations()
            # self.print_rotation_matrices(5)
            # Convert rotation matrices to rotation vectors
            self.convert_rot_mats_to_vecs()
            # self.print_rotation_vectors(5)

            k = 40
            self.clusterer: Clusterer = Clusterer(k, self.parameters, self.rotation_vecs,
                                                  self.protein_1, self.protein_2)
            self.clusterer.cluster()

            # self.protein_1.print_chain()
            # self.protein_2.print_chain()

            # self.print_residues_superimposed_results()
            # self.print_chains_superimposed_result()

        else:
            print("Unable to compare sequences.")
            return False

        return True

    def check_chain_compatibility(self):
        """
        Checks whether the residues sequence of Protein 1 and Protein 2 are similar
        :return: The percentage of similarity of the 2 protein residue sequences (Between 0 and 1)
        """
        residues_1 = [res.name for res in self.protein_1.residue_span]
        residues_2 = [res.name for res in self.protein_2.residue_span]
        sm = difflib.SequenceMatcher(None, residues_1, residues_2)
        # print(f"SM = {sm.get_matching_blocks()}")
        # print(f"Ratio = {sm.ratio()}")
        return sm.ratio()

    def check_chain_lengths(self):
        """
        Check if number of residues of Protein 1 and Protein 2 are the same.
        :return:
        """
        if self.protein_1.chain_atoms.shape[0] != self.protein_2.chain_atoms.shape[0]:
            return False
        return True

    def superimpose_chains(self):
        """
        Superimposes the entire chain of the 2 Protein structures using the backbone atoms.
        :return: A gemmi.SupResult object (Only 1, not a list)
        """
        ptype = self.protein_1.residue_span.check_polymer_type()
        self.chain_superimpose_result: gemmi.SupResult = gemmi.calculate_superposition(self.protein_1.residue_span,
                                                                                       self.protein_2.residue_span,
                                                                                       ptype, gemmi.SupSelect.MainChain)

    def sliding_window_superimpose_residues(self):
        """
        Slides a window on the backbone atoms to superimpose the residues. Produces 2 new variables:
        1. slide_window_superimpose_results = Array of superimposition information of the atoms of 2 proteins
        2. slide_window_residue_indices =   Tuple of the start and end index of one of the Protein residues array which
                                            is the span of residues after sliding a window.
        This is under the assumption that the total number of Residues of the 2 Proteins are the same.
        :return:
        """
        residues_1 = self.protein_1.residue_span
        residues_length = len(residues_1)
        backbone_1 = self.protein_1.chain_atoms.flatten()
        backbone_2 = self.protein_2.chain_atoms.flatten()
        window_size = int(self.parameters["window"])
        window_mid_index = start_index = (window_size - 1) // 2
        atom_type = self.parameters["atoms"]
        atoms_per_window = (window_size * 3) if atom_type == "backbone" else window_size

        if residues_length <= window_size:
            try:
                pos_1 = [a.pos for a in backbone_1]
                pos_2 = [a.pos for a in backbone_2]
                self.slide_window_superimpose_results.append(gemmi.superpose_positions(pos_2, pos_1))
                residues_mid_index = (residues_length - 1) // 2
                indices = (residues_mid_index, residues_mid_index + 1)
                self.protein_1.slide_window_residues_indices = indices
                self.protein_2.slide_window_residues_indices = indices
                self.get_utilised_residues(indices)
            except Exception as e:
                print(e)
            return

        try:
            final_index = residues_length + 1 - (window_size - window_mid_index)
            indices = (start_index, final_index)
            self.protein_1.slide_window_residues_indices = indices
            self.protein_2.slide_window_residues_indices = indices
            self.get_utilised_residues(indices)
            # For each residue in the 2 proteins, superimpose the backbone atoms
            for i in range(residues_length - window_size + 1):
                atoms_start_index = i * 3 if atom_type == "backbone" else i
                atoms_end_index = atoms_start_index + atoms_per_window
                # Get the x, y, and z coordinates of the atoms of the specific residue
                atoms_1 = backbone_1[atoms_start_index:atoms_end_index]
                atoms_2 = backbone_2[atoms_start_index:atoms_end_index]
                pos_1 = [a.pos for a in atoms_1]
                pos_2 = [a.pos for a in atoms_2]
                # Superimpose and append to list
                self.slide_window_superimpose_results.append(gemmi.superpose_positions(pos_2, pos_1))
        except Exception as e:
            print(e)

    def get_utilised_residues(self, indices):
        chain_1 = self.protein_1.residue_span
        chain_2 = self.protein_2.residue_span
        for i in range(indices[0], indices[1]):
            self.protein_1.slide_window_residues.append(chain_1[i])
            self.protein_2.slide_window_residues.append(chain_2[i])

    def get_transformations(self):
        """
        residues_superimpose_results is a list of gemmi.SupResult objects containing information of the superimposition
        between each residue of the Proteins. This function is used to extract the numerical data of the objects for
        KMeans clustering. Specifically, the rotation matrix.
        :return:
        """
        self.rotation_mats = np.empty(shape=[len(self.slide_window_superimpose_results), 3, 3])
        self.translation_vecs = np.empty(shape=[len(self.slide_window_superimpose_results), 3])
        for i in range(len(self.slide_window_superimpose_results)):
            rot_mat: gemmi.Mat33 = self.slide_window_superimpose_results[i].transform.mat
            trans_vec: gemmi.Vec3 = self.slide_window_superimpose_results[i].transform.vec
            self.rotation_mats[i] = rot_mat.tolist()
            self.translation_vecs[i] = trans_vec.tolist()

    def convert_rot_mats_to_vecs(self):
        self.rotation_vecs = Rotation.from_matrix(self.rotation_mats).as_rotvec(degrees=True)
        self.angles = np.array([np.linalg.norm(i) for i in self.rotation_vecs])
        self.unit_vectors = self.rotation_vecs / np.linalg.norm(self.rotation_vecs)

    # def determine_screw_axis(self):
    #     for i in range(len(self.rotation_vecs.shape[0])):
    #
    #     return

    def determine_hinge_axis(self):
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
        print(f"slide_window_residue_indices = {self.slide_window_residue_indices}")

    def print_rotation_matrices(self, n=None):
        if n is None or n > self.rotation_mats.shape[0]:
            n = self.rotation_mats.shape[0]
        print(f"rotation_mats shape = {self.rotation_mats.shape}")
        print(f"rotation_mats[0:{n}] = {self.rotation_mats[0:n]}")

    def print_rotation_vectors(self, n=None):
        if n is None or n > self.rotation_vecs.shape[0]:
            n = self.rotation_vecs.shape[0]
            for i in range(n):
                print(f"[{self.rotation_vecs[i][0]}, {self.rotation_vecs[i][1]}, {self.rotation_vecs[i][2]}],")
        print(f"rotation_vecs shape = {self.rotation_vecs.shape}")
        print(f"rotation_vecs[0:{n}] = {self.rotation_vecs[0:n]}")

    # def print_angles(self, n=None):
    #     if n is None or n > self.angles.shape[0]:
    #         n = self.angles.shape[0]
    #     print(f"angles shape = {self.angles.shape}")
    #     print(f"angles[0:{n} = {self.angles[0:n]}")

    # def superimpose_chains(self):
    #     """
    #     Superimposes each residue of the 2 Protein structures using backbone atoms.
    #     :return: List of gemmi.SupResult objects containing superimposition information
    #     """
    #     # Get the backbone atoms
    #     backbone_1 = self.protein_1.chain_atoms
    #     backbone_2 = self.protein_2.chain_atoms
    #     self.chain_superimpose_result = []
    #     try:
    #         # For each residue in the 2 proteins, superimpose the backbone atoms
    #         for i in range(backbone_1.shape[0]):
    #             # Get the x, y, and z coordinates of the backbone atoms (N, CA, C) of the specific residue
    #             pos_1 = [a.pos for a in backbone_1[i][:]]
    #             pos_2 = [a.pos for a in backbone_2[i][:]]
    #             # Superimpose and append
    #             self.chain_superimpose_result.append(gemmi.superpose_positions(pos_1, pos_2))
    #     except Exception as e:
    #         print(e)
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

