from Operations import Engine

command_file = "data/adenylate_command.txt"
data_path = "data/pdb/"
protein_structure_array = ["N", "CA", "C", "O", "S", "N", "CA", "C", "O", "H", "H"]
protein_structure_string = "NCACOSNCACOHH"


def read_command_file(file_path: str):
    new_dict = {}
    try:
        fr = open(file_path, "r")
        lines = fr.readlines()
        for line in lines:
            if not ("#" in line):
                line = line.replace("\n", "")
                tokens = line.split("=")
                new_dict[tokens[0]] = tokens[1]
        fr.close()
    except Exception as e:
        print(e)
    return new_dict


def array_sliding_window(protein_structure, window_size: int):
    """
        Possible to use yield to resume sliding window if operations are required to be done in the middle of iteration
    """
    if len(protein_structure) <= window_size:
        return protein_structure

    for i in range(len(protein_structure) - window_size + 1):
        print(protein_structure[i:i+window_size])





def main():
    command_dict = read_command_file(command_file)
    # print(command_dict)
    # array_sliding_window(protein_structure_array, int(command_dict['window']))

    file_path_1 = data_path + command_dict["filename1"]
    file_path_2 = data_path + command_dict["filename2"]

    engine = Engine(file_path_1, file_path_2, command_dict)
    # print(engine.protein_1)


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
