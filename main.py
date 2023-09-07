import FileMngr
from Operations import Engine

pdb_path = "data/pdb/"


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


if __name__ == '__main__':
    main()
