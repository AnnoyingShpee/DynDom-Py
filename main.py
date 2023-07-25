import gemmi
import Bio

NCLUSTER: int = 0
ITER: int = 0
DOMIN: int = 0
WINDLEN: int = 0
NATOT1: int = 0
NATOT2: int = 0
IAT: int = 0
IAAT: int = 0
I: int = 0
J: int = 0
NRES: int = 0
NBB: int = 0
NRES1: int = 0
NRES2: int = 0
NBB1: int = 0
NBB2: int = 0


def read_pdb_file(file):
    structure = gemmi.read_pdb(filename=file)
    return structure


def read_file(file):
    structure = gemmi.read_structure(file)
    return structure


def decode_structure(protein_structure: gemmi.Structure):
    decoded_structure = protein_structure.make_pdb_headers()
    return decoded_structure


if __name__ == '__main__':
    file_name = "pdb_files/1adg.pdb"
    # print(read_pdb_file("pdb_files/1adg.pdb"))
    structure = read_file(file_name)
    print(decode_structure(structure))
