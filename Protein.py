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
        self.structure: gemmi.Structure = gemmi.read_structure(file_path, format=gemmi.CoorFormat.Pdb)
        # self.structure.setup_entities()
        # self.structure.assign_label_seq_id()
        self.id: str = self.structure.name  # ID of the protein structure
        self.chain_param: str = chain  # The chain specified in parameter input for the program
        self.atom_type: str = atom_type
        # There is usually only one model in the structure
        self.chain: gemmi.Chain = self.structure[0][self.chain_param]
        self.chain_residues: gemmi.ResidueSpan = self.chain.get_polymer()
        self.chain_atoms = np.array(self.get_backbone_atoms())

        #########################################################
        # BioPython method
        #########################################################
        # self.id: str = file_path.replace(".pdb", "")
        # self.structure: Structure = parser.get_structure(self.id, file_path)
        # self.atom_coordinates = [atom.get_vector() for atom in self.structure.get_atoms()]
        # self.total_models, self.total_chains, self.total_residues, self.total_atoms = self.perform_protein_count()
        #########################################################

    # def __getitem__(self, item):
    #     return self.chain_atoms[]

    def get_backbone_atoms(self):
        """
        Gets the backbone atoms of each residue [(N, CA, C) or (CA only)]
        :return: A 2D array of residue atoms
        """
        atoms = []
        if self.atom_type == "backbone":
            for res in self.chain_residues:
                atoms.append([res.sole_atom(a) for a in backbone_atoms])
        elif self.atom_type == "ca":
            for res in self.chain_residues:
                atoms.append(res.sole_atom("CA"))
        elif self.atom_type == "n":
            for res in self.chain_residues:
                atoms.append(res.sole_atom("N"))
        elif self.atom_type == "c":
            for res in self.chain_residues:
                atoms.append(res.sole_atom("C"))
        return atoms

    def print_chain(self):
        print(f"{self.id}({self.chain_param}) - {self.chain_atoms.shape}")
        print(f"{self.id}({self.chain_param}) - {self.chain_atoms}")

    # def perform_protein_count(self):
    #     ##############################################
    #     # For manual counting for Gemmi
    #     ##############################################
    #     models = 0
    #     chains = 0
    #     residues = 0
    #     atoms = 0
    #     for model in self.structure:
    #         models += 1
    #         for chain in model:
    #             chains += 1
    #             for residue in chain:
    #                 residues += 1
    #                 for atom in residue:
    #                     atoms += 1
    #     ##############################################
    #
    #     ##############################################
    #     # For BioPython
    #     # https://stackoverflow.com/questions/393053/length-of-generator-output
    #     ##############################################
    #     models = sum(1 for _ in self.structure.get_models())
    #     chains = sum(1 for _ in self.structure.get_chains())
    #     residues = sum(1 for _ in self.structure.get_residues())
    #     atoms = sum(1 for _ in self.structure.get_atoms())
    #     ##############################################
    #     return models, chains, residues, atoms
