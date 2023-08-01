import sys
import gemmi
from Bio import SeqIO  # For working with sequence data (e.g., DNA, RNA, protein sequences)
from Bio.Seq import Seq  # For working with biological sequences
from Bio.SeqUtils import GC  # For calculating GC content
from Bio import AlignIO  # For handling sequence alignments
from Bio.PDB import PDBParser, PDBIO, PDBList  # For working with protein structures (PDB files)
from Bio.PDB.vectors import calc_dihedral
from Bio.PDB.Chain import Chain
from Bio.PDB.internal_coords import *
from Bio.PDB.PICIO import write_PIC, read_PIC, read_PIC_seq
from Bio.PDB.ic_rebuild import write_PDB, IC_duplicate, structure_rebuild_test
from Bio.PDB.SCADIO import write_SCAD
from Bio.SeqRecord import SeqRecord

# file = "pdb/1lfg.cif"
file = "data/pdb/1lfg.pdb"

structure: gemmi.Structure = gemmi.read_structure(file)
unit_cell: gemmi.Structure.cell = structure.cell
space_group: gemmi.Structure.spacegroup_hm = structure.spacegroup_hm
print(f"Structure =", structure)
print(f"Unit Cell =", unit_cell)
print(f"Space Group =", space_group)

# Iterate over atoms in the structure
for model in structure:  # gemmi.Model
    # print(model)  # Prints <gemmi.Model 1 with 3 chain(s)>
    for chain in model:  # gemmi.Chain
        # print(chain)  # Prints 3 <gemmi.Chain * with ** res>
        for residue in chain:  # gemmi.Residue
            # print(residue)
            # print(type(residue))
            for atom in residue:  # gemmi.Atom
                # atomic: gemmi.Atom = atom
                # print(atomic.pos)  # x, y, z coordinates as gemmi.Position
                # print(atomic.name)  # Backbone atoms as string
                # print(atomic.element)  # Atoms as gemmi.Element objects
                print(atom)
