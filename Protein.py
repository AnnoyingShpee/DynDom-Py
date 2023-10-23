import sys
import gemmi
import numpy as np
"""
Gemmi follows a hierarchy:
Structure -> Model -> Chain -> Residue -> Atom
"""

backbone_atoms = ["N", "CA", "C"]


class Protein:
    def __init__(self, file_path: str, chain: str = "A", atom_type: str = "backbone"):
        self.file_path = file_path
        self.unchanged_structure = gemmi.read_structure(file_path, format=gemmi.CoorFormat.Pdb)
        self.structure: gemmi.Structure = gemmi.read_structure(file_path, format=gemmi.CoorFormat.Pdb)
        self.id: str = self.structure.name  # ID of the protein structure
        self.chain_param: str = chain  # The chain specified in parameter input for the program
        self.atom_type: str = atom_type
        # There is usually only one model in the structure
        self.chain: gemmi.Chain = self.structure[0][self.chain_param]
        # print(f"Chain = {self.chain}")
        self.residue_span: gemmi.ResidueSpan = self.chain.get_polymer()
        # print(f"Residues = {self.chain_residues}")
        self.chain_atoms = None
        self.slide_window_residues_indices = None
        self.slide_window_residues = []
        self.unchanged_slide_window_residues = []
        self.get_backbone_atoms()

    def get_backbone_atoms(self):
        """
        Gets the backbone atoms of each residue [(N, CA, C) or (CA only)]
        :return: A 2D array of residue atoms
        """
        atoms = []
        if self.atom_type == "backbone":
            for res in self.residue_span:
                atoms.append([res.sole_atom(a) for a in backbone_atoms])
        elif self.atom_type == "ca":
            for res in self.residue_span:
                atoms.append(res.sole_atom("CA"))
        elif self.atom_type == "n":
            for res in self.residue_span:
                atoms.append(res.sole_atom("N"))
        elif self.atom_type == "c":
            for res in self.residue_span:
                atoms.append(res.sole_atom("C"))
        self.chain_atoms = np.array(atoms)

    def print_chain(self):
        print(f"{self.id}({self.chain_param}) - {self.chain_atoms.shape}")
        print(f"{self.id}({self.chain_param}) - {self.chain_atoms}")

    def print_slide_window_residues(self, n=None):
        if n is None or n > len(self.slide_window_residues):
            n = len(self.slide_window_residues)
        print(f"slide_window_residue_1 shape = {len(self.slide_window_residues)}")
        print(f"slide_window_residues_1[0:{n}] = {self.slide_window_residues[0:n]}")

    def recreate_structure(self):
        """
        To reconstruct the protein to its original structure after changes to it.
        :return:
        """
        self.slide_window_residues = []
        self.structure: gemmi.Structure = gemmi.read_structure(self.file_path, format=gemmi.CoorFormat.Pdb)
        self.chain: gemmi.Chain = self.structure[0][self.chain_param]
        self.residue_span: gemmi.ResidueSpan = self.chain.get_polymer()
        # print(f"Indices = {self.slide_window_residues_indices}")
        for i in range(self.slide_window_residues_indices[0], self.slide_window_residues_indices[1]):
            self.slide_window_residues.append(self.chain[i])



