import numpy as np
import scipy as sp
from protein.Protein import Protein
from Bio import SeqIO  # For working with sequence data (e.g., DNA, RNA, protein sequences)
from Bio.Seq import Seq  # For working with biological sequences
from Bio.SeqUtils import GC  # For calculating GC content
from Bio import AlignIO  # For handling sequence alignments
from Bio.PDB import PDBParser, PDBIO, PDBList, Superimposer  # For working with protein structures (PDB files)
from Bio.PDB.vectors import calc_dihedral
from Bio.PDB.Chain import Chain
from Bio.PDB.internal_coords import *
from Bio.PDB.PICIO import write_PIC, read_PIC, read_PIC_seq
from Bio.PDB.ic_rebuild import write_PDB, IC_duplicate, structure_rebuild_test
from Bio.PDB.SCADIO import write_SCAD
from Bio.SeqRecord import SeqRecord

import math



class Engine:
    def __init__(self, first_pdb: str, second_pdb: str,
                 params=None):
        if params is None:
            params = {"chain1id": "A",
                      "chain2id": "A",
                      "window": "5",
                      "domain": "20",
                      "ratio": "1.0"}
        self.protein_1: Protein = Protein(first_pdb)
        self.protein_2: Protein = Protein(second_pdb)
        self.parameters = params
        self.min_coordinates = []
        self.max_coordinates = []
        self.unit_vectors = np.array([])
        self.rotation_vectors = np.array([])
        self.angles = np.array([])
        self.lines = np.array([])


    def fun_rosenbrock(self, x):
        return np.array([10 * (x[1] - x[0] ** 2), (1 - x[0])])

    def least_squares_best_fit(self):
        array_1 = self.protein_1.atom_coordinates
        array_2 = self.protein_2.atom_coordinates
        centroid = np.mean(array_1, axis=0)
        cov_matrix = np.cov(array_1)
        eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
        normal_vector = eigenvectors[:, np.argmin(eigenvalues)]
        d = -np.dot(normal_vector, centroid)
        plane_equation = np.append(normal_vector, d)
        errors = np.dot(array_1, normal_vector) + d
        total_error = np.sum(errors ** 2)

        print("Best-fit plane equation:", plane_equation)
        print("Total least squares error:", total_error)
        return

    def superimpose_proteins(self):
        sup = Superimposer()


    def get_rotation_vectors(self):
        return

    def write_output_file(self):
        return

    def write_pml_file(self):
        return


