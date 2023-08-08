import sys
import gemmi
import numpy as np
import scipy as sp
# from Bio import SeqIO  # For working with sequence data (e.g., DNA, RNA, protein sequences)
# from Bio.Seq import Seq  # For working with biological sequences
# from Bio.SeqUtils import GC  # For calculating GC content
# from Bio import AlignIO  # For handling sequence alignments
# from Bio.PDB import PDBParser, PDBIO, PDBList, Structure, Model, Chain, Residue, Atom
# from Bio.PDB.vectors import calc_dihedral

"""
Gemmi follows a hierarchy:
Structure -> Model -> Chain -> Residue -> Atom


"""

"""
The Structure.cell stores the cell parameters (a, b, c, alpha, beta, gamma) 
and other properties of the cell precalculated for efficiency 
(orthogonalization and fractionalization transformations, the volume, parameters of the reciprocal unit cell).
"""

# parser = PDBParser(PERMISSIVE=1)
backbone_atoms = ["N", "CA", "C"]


class Protein:
    def __init__(self, file_path: str, chain: str = "A"):
        self.structure: gemmi.Structure = gemmi.read_structure(file_path, format=gemmi.CoorFormat.Pdb)
        # self.structure.setup_entities()
        # self.structure.assign_label_seq_id()
        self.id: str = self.structure.name  # ID of the protein structure
        self.chain_param: str = chain  # The chain specified in parameter input for the program
        # There is usually only one model in the structure
        self.chain: gemmi.Chain = self.structure[0][self.chain_param]
        self.chain_residues: gemmi.ResidueSpan = self.chain.get_polymer()
        self.chain_backbone_atoms = np.array(self.get_backbone_atoms())

        #########################################################
        # BioPython method
        #########################################################
        # self.id: str = file_path.replace(".pdb", "")
        # self.structure: Structure = parser.get_structure(self.id, file_path)
        # self.atom_coordinates = [atom.get_vector() for atom in self.structure.get_atoms()]
        # self.total_models, self.total_chains, self.total_residues, self.total_atoms = self.perform_protein_count()
        #########################################################

    def get_backbone_atoms(self):
        # atoms = np.array([])
        atoms = []
        for res in self.chain_residues:
            atoms.append([res.sole_atom(a) for a in backbone_atoms])
        return atoms

    def print_chain(self):
        print(f"{self.id}({self.chain_param}) - {self.chain_backbone_atoms.shape}")
        print(f"{self.id}({self.chain_param}) - {self.chain_backbone_atoms}")

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

    def write_to_pdb(self, path: str):
        self.structure.write_pdb(str)



# file = "data/cif/1lfg.cif"
# file = "data/pdb/2cts.pdb"
# test = Protein(file)
# print(type(test.atoms))
# for atom in test.atoms:
#     atomic: Atom.Atom = atom
#     print(atomic.coord)
# for coord in test.coordinates:
#     print(coord)
# gemmi.calculate_superposition()
# gemmi.calculate_angle()
# print(test.atom_chain)
