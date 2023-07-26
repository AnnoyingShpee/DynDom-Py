import gemmi
from Bio import SeqIO  # For working with sequence data (e.g., DNA, RNA, protein sequences)
from Bio.Seq import Seq  # For working with biological sequences
from Bio.SeqUtils import GC  # For calculating GC content
from Bio import AlignIO  # For handling sequence alignments
from Bio.PDB import PDBParser  # For working with protein structures (PDB files)
from Bio.PDB import PDBParser, PDBIO, PDBList
from Bio.PDB.vectors import calc_dihedral

pdbl = PDBList()
# pdbl.update_pdb()



def main():
    running = True
    while running:
        print("To exit application, type 'exit'.")
        structure = None
        try:
            first_input = input("Input first protein code...")
            if first_input == "exit":
                running = False
                break
            if ".pdb" in first_input:
                file_path = f"pdb_files/{first_input}"
                structure = gemmi.read_structure(file_path, format=gemmi.CoorFormat.Pdb)
            elif ".cif" in first_input:
                file_path = f"pdb_files/{first_input}"
                structure = gemmi.read_structure(file_path, format=gemmi.CoorFormat.Mmcif)
            else:
                pdb_filename = pdbl.retrieve_pdb_file(first_input)
                parser = PDBParser()
                structure = parser.get_structure(first_input, pdb_filename)
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            print(atom)
        except Exception as e:
            print(e)


if __name__ == '__main__':
    main()
