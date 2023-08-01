import sys
import gemmi
from Bio import SeqIO  # For working with sequence data (e.g., DNA, RNA, protein sequences)
from Bio.Seq import Seq  # For working with biological sequences
from Bio.SeqUtils import GC  # For calculating GC content
from Bio import AlignIO  # For handling sequence alignments
from Bio.PDB import PDBParser, PDBIO, PDBList, Structure, Model, Chain, Residue, Atom
from Bio.PDB.vectors import calc_dihedral
import numpy as np
import scipy as sp


"""
Gemmi follows a hierarchy:
Structure -> Model -> Chain -> Residue -> Atom
"""

"""
The Structure.cell stores the cell parameters (a, b, c, alpha, beta, gamma) 
and other properties of the cell precalculated for efficiency 
(orthogonalization and fractionalization transformations, the volume, parameters of the reciprocal unit cell).
"""

parser = PDBParser(PERMISSIVE=1)


class Protein:
    def __init__(self, file_path: str):
        # self.structure: gemmi.Structure = gemmi.read_structure(file_path, format=gemmi.CoorFormat.Pdb)
        self.id: str = file_path.replace(".pdb", "")
        self.structure: Structure = parser.get_structure(self.id, file_path)
        self.atom_coordinates = [atom.get_vector() for atom in self.structure.get_atoms()]
        self.total_models, self.total_chains, self.total_residues, self.total_atoms = self.perform_protein_count()

    # def initialise_variables(self):
    #     # chain_array, coordinates = numpy.array([]), numpy.array([])
    #     for model in self.structure:
    #         for chain in model:
    #             for residue in chain:
    #                 for atom in residue:
    #                     atomic: gemmi.Atom = atom
    #                     self.atoms = np.append(self.atoms, atomic.name)
    #                     pos: gemmi.Position = atomic.pos
    #                     self.coordinates = np.append(self.coordinates, (pos.x, pos.y, pos.z))
    #                     # pos_norm: gemmi.Vec3 = pos.normalized()
    #                     # self.coordinates_norm = np.append(self.coordinates, (pos_norm.x, pos_norm.y, pos_norm.z))

    def perform_protein_count(self):
        # https://stackoverflow.com/questions/393053/length-of-generator-output
        # models = 0
        # chains = 0
        # residues = 0
        # atoms = 0
        # for model in self.structure:
        #     models += 1
        #     for chain in model:
        #         chains += 1
        #         for residue in chain:
        #             residues += 1
        #             for atom in residue:
        #                 atoms += 1
        models = sum(1 for _ in self.structure.get_models())
        chains = sum(1 for _ in self.structure.get_chains())
        residues = sum(1 for _ in self.structure.get_residues())
        atoms = sum(1 for _ in self.structure.get_atoms())
        return models, chains, residues, atoms

    def write_to_pdb(self, path: str):
        self.structure.write_pdb(str)

    # def get_cif_document(self):
    #     self.structure: gemmi.cif.Document = gemmi.cif.read_file(self.file_name)
    #     print(self.structure.as_string())
    #     # self.structure = gemmi.read_structure(self.file_name, format=gemmi.CoorFormat.Mmcif)
    #     # Returns a gemmi.cif.Block object. Can
    #     # self.structure_metadata: gemmi.cif.Block = self.structure.make_mmcif_headers()

    # def get_pdb_structure(self):
    #     self.structure = gemmi.read_structure(self.file_name, format=gemmi.CoorFormat.Pdb)
    #     # Returns a string of the whole structure's details
    #     self.structure_metadata = self.structure.make_pdb_headers()

    # def get_structure(self):
    #     self.structure = gemmi.read_structure(self.file_name, format=gemmi.CoorFormat.Detect)
    #     self.structure_metadata = self.structure
    #
    # def print_pipeline(self):
    #     if self.file_type == "pdb":
    #         self.get_pdb_structure()
    #     elif self.file_type == "cif":
    #         self.get_cif_document()

    # https://gemmi.readthedocs.io/en/latest/mol.html#mcra

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
