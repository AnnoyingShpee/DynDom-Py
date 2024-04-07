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
                param_name = tokens[0]
                param_val = tokens[1]
                if param_name == "window":
                    param_val = int(tokens[1])
                elif param_name == "domain":
                    param_val = int(tokens[1])
                elif param_name == "ratio":
                    param_val = float(tokens[1])
                temp_dict[param_name] = param_val
        fr.close()
    except Exception as e:
        print(e)
    return temp_dict


def write_rotation_vec_to_pdb(file_name: str, slide_window_residues, slide_window_indices, rotation_vectors):
    """
    Writes the rotation vectors of each residue of the slide window into a pdb file
    :param file_name:
    :param slide_window_residues:
    :param slide_window_indices:
    :param rotation_vectors:
    :return:
    """
    try:
        fw = open(f"{output_pdb_file_path}{file_name}_rot_vecs.pdb", "w")
        start_index = slide_window_indices[0]
        for r in range(len(slide_window_residues)):
            residue_name = slide_window_residues[r].name
            residue_num = str(start_index + r + 1).rjust(3, " ")
            x = str(round(rotation_vectors[r][0], 3)).rjust(8, " ")
            y = str(round(rotation_vectors[r][1], 3)).rjust(8, " ")
            z = str(round(rotation_vectors[r][2], 3)).rjust(8, " ")
            row = f"ATOM         CA  {residue_name} A {residue_num}    {x}{y}{z}\n"
            fw.write(row)
        fw.close()
    except Exception as e:
        print(e)
        return False
    return True


def write_w5_info_file(file_name: str, protein_1_name:str, protein_2_name: str, param: dict, domains: list,
                       fixed_domain_id: int):
    try:
        fw = open(f"{protein_1_name}{param['chain1id']}_{protein_2_name}{param['chain2id']}.w5_info", "w")
        fw.write("DynDom Python Version 1.0\n")
        fw.write(f"{protein_1_name}{param['chain1id']}_{protein_2_name}{param['chain2id']}.w5\n")
        fw.write(f"file name of conformer 1: {protein_1_name}.pdb\n")
        fw.write(f"chain id: {param['chain1id']}\n")
        fw.write(f"file name of conformer 2: {protein_2_name}.pdb\n")
        fw.write(f"chain id: {param['chain2id']}\n")
        fw.write(f"window length: {param['window']}\n")
        fw.write(f"minimum ratio of external to internal motion: {param['ratio']}\n")
        fw.write(f"minimum domain size: {param['domain']}\n")
        fw.write(f"atoms to use: {param['atoms']}\n")
        fw.write(f"THERE ARE {len(domains)} DOMAINS\n")
        fw.write("================================================================================\n")
        fw.write("FIXED DOMAIN\n")
        fw.write(f"DOMAIN NUMBER: \t {fixed_domain_id} (coloured blue for rasmol)\n")
        fw.write(f"RESIDUE NUMBERS: ")

    except Exception as e:
        print(e)
        return False
    return True



def write_pymol_file(pml_file_name: str, pdb_file_name: str, data):
    try:
        fw = open(f"{output_pymol_file_path}{pml_file_name}.pml", "w")
        fw.write("reinitialize\n")
        fw.write(f"load {pdb_file_name}\n")
        fw.write("bg_color white")
        fw.write("color grey")
        fw.close()
    except Exception as e:
        print(e)
        return False
    return True

