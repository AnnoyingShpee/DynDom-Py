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
                               "ratio": "1.0"}
        else:
            self.parameters = params
        # Initialise proteins
        self.protein_1: Protein = Protein(first_pdb, self.parameters["chain1id"])
        self.protein_2: Protein = Protein(second_pdb, self.parameters["chain2id"])
        # Array of gemmi.SupResult objects containing superimposition information between each residue
        self.residues_superimpose_results = np.array([])
        # A gemmi.SupResult object containing superimposition information between 2 chains
        self.chain_superimpose_result = None
        # List of arrays of gemmi.SupResult objects containing superimposition information between each residue
        self.slide_window_superimpose_results = np.array([])
        # List of (3 x 3) rotation matrices
        self.rotation_mats = np.empty(shape=(self.slide_window_superimpose_results.shape[0], 3, 3))
        # List of rotation vectors created from rotation_mats
        self.rotation_vecs = np.empty(shape=(self.rotation_mats.shape[0], 3))
        # List of angles
        self.angles = np.empty(shape=(1, self.rotation_mats.shape[0]))
        # self.slide_window_result_1, self.slide_window_result_2 = np.array([]), np.array([])

    def run(self):
        if self.check_chain_compatibility() > 0.4 and self.check_chain_lengths():
            self.superimpose_chains()
            self.sliding_window_superimpose()

            self.get_rotation_mats()
            self.convert_rot_mats_to_vecs()

            k = 40

            clusterer = Clusterer(k, self.parameters, self.rotation_vecs)
            clusterer.cluster()


            # self.slide_window_result_1, self.slide_window_result_2 = self.sliding_window_on_backbone_atoms_1d()
            # self.slide_window_result_1, self.slide_window_result_2 = self.sliding_window_on_backbone_atoms_2d()
            # self.print_slide_window_result()

            # Clusterer.calc_k_means_sklearn(self.rotation_vecs, 25)



            # FileMngr.write_pdb_file()
            # FileMngr.write_pymol_file()
        # else:
        #     print("Sequences are too different to compare")
        #     return False

        # self.protein_1.print_chain()
        # self.protein_2.print_chain()

        # self.print_residues_superimposed_results()
        # self.print_chains_superimposed_result()

        self.print_slide_window_superimpose_result(5)

        # self.print_rotation_matrices(5)
        # self.print_rotation_vectors(5)
        # self.print_angles(5)

        return True

    def check_chain_compatibility(self):
        """
        Checks whether the residues sequence of Protein 1 and Protein 2 are similar
        :return: The percentage of similarity of the 2 protein residue sequences (Between 0 and 1)
        """
        residues_1 = [res.name for res in self.protein_1.chain_residues]
        residues_2 = [res.name for res in self.protein_2.chain_residues]
        sm = difflib.SequenceMatcher(None, residues_1, residues_2)
        # print(f"SM = {sm.get_matching_blocks()}")
        # print(f"Ratio = {sm.ratio()}")
        return sm.ratio()

    def check_chain_lengths(self):
        """
        Check if number of residues of Protein 1 and Protein 2 are the same.
        :return:
        """
        if self.protein_1.chain_backbone_atoms.shape[0] != self.protein_2.chain_backbone_atoms.shape[0]:
            return False
        return True

    def superimpose_residues(self):
        """
        Superimposes each residue of the 2 Protein structures using backbone atoms.
        :return: List of gemmi.SupResult objects containing superimposition information
        """
        # Get the backbone atoms
        backbone_1 = self.protein_1.chain_backbone_atoms
        backbone_2 = self.protein_2.chain_backbone_atoms
        try:
            # For each residue in the 2 proteins, superimpose the backbone atoms
            for i in range(backbone_1.shape[0]):
                # Get the x, y, and z coordinates of the backbone atoms (N, CA, C) of the specific residue
                pos_1 = [a.pos for a in backbone_1[i][:]]
                pos_2 = [a.pos for a in backbone_2[i][:]]
                # Superimpose and append
                self.residues_superimpose_results = np.append(self.residues_superimpose_results,
                                                              gemmi.superpose_positions(pos_1, pos_2))
        except Exception as e:
            print(e)

    def superimpose_chains(self):
        """
        Superimposes the entire chain of the 2 Protein structures using the backbone atoms.
        :return: A gemmi.SupResult object (Only 1, not a list)
        """
        polymer_1 = self.protein_1.chain_residues
        polymer_2 = self.protein_2.chain_residues
        ptype = polymer_1.check_polymer_type()
        self.chain_superimpose_result: gemmi.SupResult = gemmi.calculate_superposition(polymer_1, polymer_2,
                                                                                       ptype, gemmi.SupSelect.MainChain)

    def sliding_window_superimpose(self):
        """
        Slides a window on the backbone atoms to superimpose the residues.
        This is under the assumption that the number of backbone atoms (Residues) of the 2 Proteins are the same.
        :return:
        """
        backbone_1 = self.protein_1.chain_backbone_atoms.flatten()
        backbone_2 = self.protein_2.chain_backbone_atoms.flatten()
        window_size = int(self.parameters["window"])

        if len(backbone_1) <= window_size:
            try:
                # For each residue in the 2 proteins, superimpose the backbone atoms
                for i in range(backbone_1.shape[0]):
                    # Get the x, y, and z coordinates of the backbone atoms (N, CA, C) of the specific residue
                    pos_1 = [a.pos for a in backbone_1[i][:]]
                    pos_2 = [a.pos for a in backbone_2[i][:]]
                    # Superimpose and append to array
                    self.slide_window_superimpose_results = np.append(self.slide_window_superimpose_results,
                                                                      gemmi.superpose_positions(pos_1, pos_2))
            except Exception as e:
                print(e)
            return

        try:
            for i in range(len(backbone_1) - window_size + 1):
                contents_1 = backbone_1[i:i + window_size]
                contents_2 = backbone_2[i:i + window_size]
                pos_1 = [a.pos for a in contents_1]
                pos_2 = [a.pos for a in contents_2]
                self.slide_window_superimpose_results = np.append(self.slide_window_superimpose_results,
                                                                  gemmi.superpose_positions(pos_1, pos_2))
        except Exception as e:
            print(e)

    def get_rotation_mats(self):
        """
        residues_superimpose_results is a list of gemmi.SupResult objects containing information of the superimposition
        between each residue of the Proteins. This function is used to extract the numerical data of the objects for
        KMeans clustering. Specifically, the rotation matrix.
        :return:
        """
        temp = []
        for i in range(len(self.slide_window_superimpose_results)):
            matrix: gemmi.Mat33 = self.slide_window_superimpose_results[i].transform.mat
            listed_mat = matrix.tolist()
            temp.append(listed_mat)
            # mat_1d = np.array(listed_mat).flatten()
            # temp.append(mat_1d)
        self.rotation_mats = np.array(temp)

    def convert_rot_mats_to_vecs(self):
        rot_mats = [Rotation.from_matrix(sr.transform.mat.tolist()) for sr in self.slide_window_superimpose_results]
        self.rotation_vecs = np.array([rm.as_rotvec(degrees=True) for rm in rot_mats])
        self.angles = np.array([np.linalg.norm(i) for i in self.rotation_vecs])

    def print_residues_superimposed_results(self, n=None):
        if n is None or n > len(self.residues_superimpose_results):
            n = len(self.residues_superimpose_results)
        print(f"Superimposed residues results shape = {self.residues_superimpose_results.shape}")
        for i in range(n):
            result_obj: gemmi.SupResult = self.residues_superimpose_results[i]
            print(f"{i+1} RMSD = {result_obj.rmsd}")
            print(f"    Count (Number of atoms used in each chain) = {result_obj.count}")
            print(f"    Center 1 = {result_obj.center1}")
            print(f"    Center 2 = {result_obj.center2}")
            print(f"    Translation Vector = {result_obj.transform.vec}")
            print(f"    Rotation Matrix = {result_obj.transform.mat}")

    def print_chains_superimposed_result(self):
        print(f"RMSD =                  {self.chain_superimpose_result.rmsd}")
        print(f"Count =                 {self.chain_superimpose_result.count}")
        print(f"Center 1 =              {self.chain_superimpose_result.center1}")
        print(f"Center 2 =              {self.chain_superimpose_result.center2}")
        print(f"Translation Vector =    {self.chain_superimpose_result.transform.vec}")
        print(f"Rotation Matrix =       {self.chain_superimpose_result.transform.mat}")

    # def print_slide_window_result(self, n=None):
    #     if n is None or n > self.slide_window_result_1.shape[0]:
    #         n = self.slide_window_result_1.shape[0]
    #     print(f"slide_window_result_1 shape = {self.slide_window_result_1.shape}")
    #     print(f"slide_window_result_1[0:{n}] = {self.slide_window_result_1[0:n]}")
    #     print(f"slide_window_result_2 shape = {self.slide_window_result_2.shape}")
    #     print(f"slide_window_result_2[0:{n}] = {self.slide_window_result_2[0:n]}")

    def print_slide_window_superimpose_result(self, n=None):
        if n is None or n > self.slide_window_superimpose_results.shape[0]:
            n = self.slide_window_superimpose_results.shape[0]
        print(f"slide_window_superimpose_result shape = {self.slide_window_superimpose_results.shape}")
        print(f"slide_window_superimpose_result[0:{n}] = {self.slide_window_superimpose_results[0:n]}")
        for i in range(n):
            item = self.slide_window_superimpose_results[i]
            print(f"RMSD =                  {item.rmsd}")
            print(f"Count =                 {item.count}")
            print(f"Center 1 =              {item.center1}")
            print(f"Center 2 =              {item.center2}")
            print(f"Translation Vector =    {item.transform.vec}")
            print(f"Rotation Matrix =       {item.transform.mat}")

    def print_rotation_matrices(self, n=None):
        if n is None or n > self.rotation_mats.shape[0]:
            n = self.rotation_mats.shape[0]
        print(f"rotation_mats shape = {self.rotation_mats.shape}")
        print(f"rotation_mats[0:{n}] = {self.rotation_mats[0:n]}")

    def print_rotation_vectors(self, n=None):
        if n is None or n > self.rotation_vecs.shape[0]:
            n = self.rotation_vecs.shape[0]
        print(f"rotation_vecs shape = {self.rotation_vecs.shape}")
        print(f"rotation_vecs[0:{n}] = {self.rotation_vecs[0:n]}")

    def print_angles(self, n=None):
        if n is None or n > self.angles.shape[0]:
            n = self.angles.shape[0]
        print(f"angles shape = {self.angles.shape}")
        print(f"angles[0:{n} = {self.angles[0:n]}")







    # Creates a 1D array of the coordinates
    # def sliding_window_on_backbone_atoms_1d(self):
    #     """
    #     Slides a window on the backbone atom chains to obtain the x, y, and z coordinates.
    #     Sliding window has a step of 1 when sliding across the chains.
    #     This is under the assumption that the number of residues of the 2 proteins are the same.
    #     :return: 1D array of the coordinates of the windowed backbone atoms.
    #     """
    #     # chain_backbone_atoms is a 2D array, hence, it needs to be converted into a 1D array
    #     backbone_1 = self.protein_1.chain_backbone_atoms.flatten()
    #     backbone_2 = self.protein_2.chain_backbone_atoms.flatten()
    #     temp_1 = np.array([])
    #     temp_2 = np.array([])
    #     # Get the window size
    #     window_size = int(self.parameters["window"])
    #     # If the total number of backbone atoms are less than the window size, no need to slide a window.
    #     if len(backbone_1) <= window_size:
    #         return np.array(backbone_1), np.array(backbone_2)
    #     try:
    #         for i in range(len(backbone_1) - window_size + 1):
    #             temp_1 = np.append(temp_1, backbone_1[i:i + window_size])
    #             temp_2 = np.append(temp_2, backbone_2[i:i + window_size])
    #     except Exception as e:
    #         print(e)
    #     return temp_1, temp_2

    # def sliding_window_on_backbone_atoms_2d(self):
    #     """
    #     lides a window on the backbone atom chains to obtain the x, y, and z coordinates.
    #     Sliding window has a step of 1 when sliding across the chains.
    #     This is under the assumption that the number of residues of the 2 proteins are the same.
    #     :return: 2D array of the coordinates [N, window_size]
    #     """
    #     # chain_backbone_atoms is a 2D array, hence, it needs to be converted into a 1D array
    #     backbone_1 = self.protein_1.chain_backbone_atoms.flatten()
    #     backbone_2 = self.protein_2.chain_backbone_atoms.flatten()
    #     # Get the window size
    #     window_size = int(self.parameters["window"])
    #     temp_1 = []
    #     temp_2 = []
    #     # If the total number of backbone atoms are less than the window size, no need to slide a window.
    #     if len(backbone_1) <= window_size:
    #         return np.array(backbone_1), np.array(backbone_2)
    #     try:
    #         for i in range(len(backbone_1) - window_size + 1):
    #             temp_1.append(backbone_1[i:i + window_size])
    #             temp_2.append(backbone_2[i:i + window_size])
    #     except Exception as e:
    #         print(e)
    #
    #     return np.array(temp_1), np.array(temp_2)