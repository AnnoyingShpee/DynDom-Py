import sys
input_command_file_path = "data/"
input_pdb_file_path = "data/pdb/"
output_pymol_file_path = "output/pml/"
output_pdb_file_path = "output/pdb/"


def read_command_file(file_name: str):
    temp_dict = {}
    try:
        fr = open(f"{input_command_file_path}{file_name}", "r")
        lines = fr.readlines()
        for line in lines:
            if not ("#" in line):
                line = line.replace("\n", "")
                line = line.replace(" ", "")
                tokens = line.split("=")
                temp_dict[tokens[0]] = tokens[1]
        fr.close()
    except Exception as e:
        print(e)
    return temp_dict


def write_pdb_file(file_name: str, data):
    try:
        fw = open(f"{output_pdb_file_path}{file_name}", "w")
    except Exception as e:
        print(e)
        return False
    return True


def write_pymol_file(file_name: str, data):
    try:
        fw = open(f"{output_pymol_file_path}{file_name}", "w")

    except Exception as e:
        print(e)
        return False
    return True
