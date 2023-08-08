import sys
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

import matplotlib.pyplot as plt


# # Seq object different from Python string
# test_seq = Seq("AGTACACTGGT")
# test_string = "AGTACACTGGT"
#
# print(test_seq)
# print(test_seq.complement())  # TCATGTGACCA
# print(test_seq.reverse_complement())  # ACCAGTGTACT
#
#
# pdb_file = "pdb/1lfg.pdb"
# pdb_parser = PDBParser()
# structure = pdb_parser.get_structure("A Protein", pdb_file)
# print(structure)
# header = structure.header
# print(header)
# name = structure.header["name"]
# resolution = structure.header["resolution"]
# keywords = structure.header["keywords"]
# print(keywords)

# fasta_file = "fasta_files/ls_orchid.fasta"
# for seq_record in SeqIO.parse(fasta_file, "fasta"):
#     print(seq_record.id)
#     print(repr(seq_record.seq))
#     print(len(seq_record))

# gbk_file = "gbk_files/ls_orchid.gbk"
# for seq_record in SeqIO.parse(gbk_file, "genbank"):
#     print(seq_record.id)
#     print(repr(seq_record.seq))
#     print(len(seq_record))


def main():
    # Download PDB file (if not already present)
    pdb_code = "1CRN"  # Replace with your PDB code
    pdbl = PDBList()
    # pdbl = PDBList(pdb="data/pdb")
    # pdbl.update_pdb()
    pdb_filename = pdbl.retrieve_pdb_file(pdb_code)

    # Parse PDB file
    parser = PDBParser()
    structure = parser.get_structure(pdb_code, pdb_filename)



if __name__ == "__main__":
    main()
