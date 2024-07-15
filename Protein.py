import sys
import gemmi
import numpy as np
"""
Gemmi follows a hierarchy:
Structure -> Model -> Chain -> Residue -> Atom
"""


class Protein:
    def __init__(self, file_path: str, chain: str = "A", atom_type: str = "backbone"):
        self.file_path = file_path
        self.structure: gemmi.Structure = gemmi.read_structure(file_path, format=gemmi.CoorFormat.Pdb)
        self.id: str = self.structure.name  # ID of the protein structure
        self.chain_param: str = chain  # The chain specified in parameter input for the program
        if atom_type == "backbone":
            self.atoms_to_use = ["N", "CA", "C"]
        elif atom_type == "ca":
            self.atoms_to_use = ["CA"]
        # There is usually only one model in the structure
        self.chain: gemmi.Chain = self.structure[0][self.chain_param]
        self.residue_span: gemmi.ResidueSpan = self.chain.get_polymer()
        self.utilised_residues_indices = []
        self.chain_atoms = None
        self.slide_window_residues_indices = None

    def get_backbone_atoms(self):
        """
        Gets the backbone atoms of each residue [(N, CA, C) or (CA only)]
        :return: A 2D array of residue atoms
        """
        atoms = []
        for i in self.utilised_residues_indices:
            atoms.append([self.residue_span[i].sole_atom(a) for a in self.atoms_to_use])
        self.chain_atoms = np.asarray(atoms)

    def get_structure(self):
        return gemmi.read_structure(self.file_path, format=gemmi.CoorFormat.Pdb)

    def get_model(self):
        structure = gemmi.read_structure(self.file_path, format=gemmi.CoorFormat.Pdb)
        return structure[0]

    def get_chain(self):
        structure = gemmi.read_structure(self.file_path, format=gemmi.CoorFormat.Pdb)
        return structure[0][self.chain_param]

    def get_polymer(self):
        structure = gemmi.read_structure(self.file_path, format=gemmi.CoorFormat.Pdb)
        return structure[0][self.chain_param].get_polymer()

    def get_utilised_residues(self):
        chain: gemmi.ResidueSpan = self.get_polymer()
        slide_window_chain = gemmi.Chain(self.chain_param)
        for i in range(len(self.utilised_residues_indices)):
            index = self.utilised_residues_indices[i]
            slide_window_chain.add_residue(chain[index])
        return slide_window_chain.get_polymer()

    def get_slide_window_residues(self):
        """
        Get the residues from only the sliding window result.
        :return:
        """
        chain: gemmi.ResidueSpan = self.get_polymer()
        slide_window_chain = gemmi.Chain(self.chain_param)
        for i in range(self.slide_window_residues_indices[0], self.slide_window_residues_indices[1]):
            index = self.utilised_residues_indices[i]
            slide_window_chain.add_residue(chain[index])
        return slide_window_chain.get_polymer()

    def print_chain(self):
        print(f"{self.id}({self.chain_param}) - {self.chain_atoms.shape}")
        print(f"{self.id}({self.chain_param}) - {self.chain_atoms}")

