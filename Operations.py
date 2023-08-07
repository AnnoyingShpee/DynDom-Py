import numpy as np
import scipy as sp
import gemmi
from protein.Protein import Protein
import math

# rotation_matrix_x = np.array([[1,     0,              0],
#                               [0,     math.cos(),     -math.sin()],
#                               [0,     math.sin(),     math.cos()]])
#
# rotation_matrix_y = np.array([[math.cos(),    0,      math.sin()],
#                               [0,             1,      0],
#                               [-math.sin(),   0,      math.cos()]])
#
# rotation_matrix_z = np.array([[math.cos(),    -math.sin(),        0],
#                               [math.sin(),    math.cos(),         0],
#                               [0,             0,                  1]])


class Engine:
    def __init__(self, first_pdb: str, second_pdb: str, params=None):
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
        self.superimpose_results = np.array([])
        self.slide_window_result_1, self.slide_window_result_2 = np.array([]), np.array([])
        self.min_coordinates = []
        self.max_coordinates = []
        self.unit_vectors = np.array([])
        self.rotation_vectors = np.array([])
        self.angles = np.array([])
        self.lines = np.array([])

    def run(self):
        self.superimpose_residues()
        self.sliding_window_on_backbone_atoms()
        self.superimpose_polymers()
        # self.protein_1.print_chain()
        # self.protein_2.print_chain()
        # print(f"Superimposed results shape = {self.superimpose_results.shape}")
        # for i in range(len(self.superimpose_results)):
        #     result_obj: gemmi.SupResult = self.superimpose_results[i]
        #     print(f"{i+1} RMSD = {result_obj.rmsd}")
        #     print(f"    Count (Number of atoms used in each chain) = {result_obj.count}")
        #     print(f"    Center 1 = {result_obj.center1}")
        #     print(f"    Center 2 = {result_obj.center2}")
        #     print(f"    Translation Vector = {result_obj.transform.vec}")
        #     print(f"    Rotation Matrix = {result_obj.transform.mat}")
        # print(f"slide_window_result_1 shape = {self.slide_window_result_1.shape}")
        # print(f"slide_window_result_1[0:6] = {self.slide_window_result_1[0:6]}")
        # self.center_1, self.center_2 = self.superimpose_chains()
        # self.slide_window_result_1, self.slide_window_result_2 = self.sliding_window_on_backbone_atoms()
        return

    def superimpose_residues(self):
        backbone_1 = self.protein_1.chain_backbone_atoms
        backbone_2 = self.protein_2.chain_backbone_atoms
        if backbone_1.shape[0] != backbone_2.shape[0]:
            print("Chains do not have the same number of residues")
            return False
        try:
            for i in range(backbone_1.shape[0]):
                pos_1 = [a.pos for a in backbone_1[i][:]]
                pos_2 = [a.pos for a in backbone_2[i][:]]
                self.superimpose_results = np.append(self.superimpose_results, gemmi.superpose_positions(pos_1, pos_2))
        except Exception as e:
            print(e)
            return False
        return True

    def superimpose_polymers(self):
        polymer_1 = self.protein_1.chain_residues
        polymer_2 = self.protein_2.chain_residues
        ptype = polymer_1.check_polymer_type()
        # print(ptype)
        sup: gemmi.SupResult = gemmi.calculate_superposition(polymer_1, polymer_2, ptype, gemmi.SupSelect.MainChain)
        print(f"RMSD =                  {sup.rmsd}")
        print(f"Count =                 {sup.count}")
        print(f"Center 1 =              {sup.center1}")
        print(f"Center 2 =              {sup.center2}")
        print(f"Translation Vector =    {sup.transform.vec}")
        print(f"Rotation Matrix =       {sup.transform.mat}")


    # Creates a 1D array of the coordinates
    def sliding_window_on_backbone_atoms(self):
        backbone_1 = self.protein_1.chain_backbone_atoms.flatten()
        backbone_2 = self.protein_2.chain_backbone_atoms.flatten()
        window_size = int(self.parameters["window"])
        if len(backbone_1) != len(backbone_2):
            print("Not the same size")
            return False

        if len(backbone_1) <= window_size:
            self.slide_window_result_1, self.slide_window_result_2 = np.array(backbone_1), np.array([backbone_2])
            return True

        try:
            for i in range(len(backbone_1) - window_size + 1):
                self.slide_window_result_1 = np.append(self.slide_window_result_1, backbone_1[i:i + window_size])
                self.slide_window_result_2 = np.append(self.slide_window_result_2, backbone_2[i:i + window_size])
        except Exception as e:
            print(e)
            return False

        return True

    # # Creates a 2D array of the coordinates [N, window_size]
    # def sliding_window_on_backbone_atoms(self):
    #     backbone_1 = self.protein_1.chain_backbone_atoms.flatten()
    #     backbone_2 = self.protein_2.chain_backbone_atoms.flatten()
    #     window_size = int(self.parameters["window"])
    #     array_1 = []
    #     array_2 = []
    #     if len(backbone_1) != len(backbone_2):
    #         print("Not the same size")
    #         return False
    #
    #     if len(backbone_1) <= window_size:
    #         self.slide_window_result_1, self.slide_window_result_2 = np.array(backbone_1), np.array([backbone_2])
    #         return True
    #
    #     try:
    #         for i in range(len(backbone_1) - window_size + 1):
    #             array_1.append(backbone_1[i:i + window_size])
    #             array_2.append(backbone_2[i:i + window_size])
    #     except Exception as e:
    #         print(e)
    #         return False
    #
    #     self.slide_window_result_1, self.slide_window_result_2 = np.array(array_1), np.array([array_2])
    #     return True

    def get_rotation_vectors(self):
        return

    def write_result_file(self):
        return

    def write_pml_file(self):
        return


