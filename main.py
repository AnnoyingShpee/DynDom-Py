import FileMngr
from Operations import Engine

command_file_path = "data/adenylate_command.txt"
pdb_path = "data/pdb/"


def array_sliding_window(protein_structure, window_size: int):
    """
        Possible to use yield to resume sliding window if operations are required to be done in the middle of iteration
    """
    if len(protein_structure) <= window_size:
        return protein_structure

    for i in range(len(protein_structure) - window_size + 1):
        print(protein_structure[i:i+window_size])


def main():
    # Read command file to get parameters ( Protein PDB file names, protein chains, window size, domain size, ratio )
    param_dict = FileMngr.read_command_file("adenylate_command.txt")
    # Concatenate PDB file names with path
    pdb_path_1 = pdb_path + param_dict["filename1"]
    pdb_path_2 = pdb_path + param_dict["filename2"]
    # Initialise Engine object
    engine = Engine(pdb_path_1, pdb_path_2, param_dict)
    # Run the Engine
    engine.run()

    # running = True
    # while running:
    #     print("To exit application, type 'exit'.")
    #     structure = None
    #     try:
    #         first_input = input("Input first protein code...")
    #         if first_input == "exit":
    #             running = False
    #             break
    #         if ".pdb" in first_input:
    #             file_path = f"pdb_files/{first_input}"
    #             structure = gemmi.read_structure(file_path, format=gemmi.CoorFormat.Pdb)
    #         elif ".cif" in first_input:
    #             file_path = f"pdb_files/{first_input}"
    #             structure = gemmi.read_structure(file_path, format=gemmi.CoorFormat.Mmcif)
    #         else:
    #             pdb_filename = pdbl.retrieve_pdb_file(first_input)
    #             parser = PDBParser()
    #             structure = parser.get_structure(first_input, pdb_filename)
    #         for model in structure:
    #             for chain in model:
    #                 for residue in chain:
    #                     for atom in residue:
    #                         print(atom)
    #     except Exception as e:
    #         print(e)


if __name__ == '__main__':
    main()
