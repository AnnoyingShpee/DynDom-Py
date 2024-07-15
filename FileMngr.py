import traceback
import urllib.request
from itertools import groupby
from operator import itemgetter
from pathlib import Path

input_command_file_path = "data"
input_pdb_file_path = "data/pdb"
output_path = "output"


def read_command_file():
    temp_dict = {}
    try:
        fr = open(f"{input_command_file_path}/command.txt", "r")
        lines = fr.readlines()
        for line in lines:
            if not ("#" in line):
                line = line.replace("\n", "")
                line = line.replace(" ", "")
                tokens = line.split("=")
                param_name = tokens[0]
                param_val = tokens[1]
                if "chain" in param_name:
                    param_val = param_val.upper()
                elif "filename" in param_name:
                    param_val = param_val.lower()
                temp_dict[param_name] = param_val
        fr.close()
        make_sure_pdb_exists(temp_dict["filename1"])
        make_sure_pdb_exists(temp_dict["filename2"])
    except Exception as e:
        print(e)
    return temp_dict


def make_sure_pdb_exists(file_name):
    """
    Uses HTTPS request to make sure that the pdb file is in the directory
    :param file_name:
    :return:
    """
    file_path = f"{input_pdb_file_path}/{file_name}"
    path = Path(file_path)
    if not path.exists():
        try:
            urllib.request.urlretrieve(
                f"https://files.rcsb.org/download/{file_name}",
                file_path
            )
        except Exception as e:
            print(e)


def read_param_file():
    temp_dict = {}
    try:
        fr = open(f"{input_command_file_path}/param.txt", "r")
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
                elif param_name == "k_means_n_init":
                    param_val = int(tokens[1])
                elif param_name == "k_means_max_iter":
                    param_val = int(tokens[1])
                temp_dict[param_name] = param_val
        fr.close()
    except Exception as e:
        print(e)
    return temp_dict


def write_rotation_vec_to_pdb(folder_name: str, protein, rotation_vectors):
    """
    Writes the rotation vectors of each residue of the slide window into a pdb file
    :param folder_name:
    :param protein
    :param rotation_vectors:
    :return:
    """
    try:
        dir_path_str = f"{output_path}/{folder_name}"
        dir_path = Path(dir_path_str)
        if not dir_path.exists():
            dir_path.mkdir(parents=True)
        fw = open(f"{dir_path_str}/{protein.id}_rot_vecs.pdb", "w")
        slide_window_indices = protein.slide_window_residues_indices
        protein_polymer = protein.get_polymer()
        for i in range(rotation_vectors.shape[0]):
            index = protein.utilised_residues_indices[i+slide_window_indices[0]]
            residue_name = protein_polymer[index].name
            residue_num = protein_polymer[index].seqid.num
            x = str(round(rotation_vectors[i][0], 3)).rjust(8, " ")
            y = str(round(rotation_vectors[i][1], 3)).rjust(8, " ")
            z = str(round(rotation_vectors[i][2], 3)).rjust(8, " ")
            row = f"ATOM         CA  {residue_name} A {residue_num}    {x}{y}{z}\n"
            fw.write(row)
        fw.close()
    except Exception as e:
        print(e)
        return False
    return True


def write_final_output_pdb(folder_name, protein_1, fitted_protein_2, fitted_protein_2_chain):
    """
    :param folder_name: The name of the output folder
    :param protein_1: The Protein object of protein 1
    :param fitted_protein_2: The ResidueSpan object of protein 2 fitted to protein 1
    :param fitted_protein_2_chain: The chain of protein 2 fitted to protein 1
    :return:
    """
    try:
        dir_path_str = f"{output_path}/{folder_name}"
        dir_path = Path(dir_path_str)
        if not dir_path.exists():
            dir_path.mkdir(parents=True)
        fw = open(f"{dir_path_str}/{folder_name}.pdb", "w")
        fw.write(f"MODEL{'1'.rjust(9, ' ')}\n")
        protein_1_residues = protein_1.get_polymer()
        atom_count = 1
        subchain = protein_1.chain_param
        for r in protein_1_residues:
            res_name = r.name.rjust(3, " ")
            res_num = str(r.seqid.num).rjust(4, " ")
            for a in r:
                atom_num = str(atom_count).rjust(5, " ")
                atom_name = a.name.ljust(4, " ")
                x = str(round(a.pos.x, 3)).rjust(8, " ")
                y = str(round(a.pos.y, 3)).rjust(8, " ")
                z = str(round(a.pos.z, 3)).rjust(8, " ")
                row = f"ATOM  {atom_num} {atom_name} {res_name} {subchain}{res_num}    {x}{y}{z}\n"
                fw.write(row)
                atom_count += 1
        fw.write("ENDMDL\n")

        fw.write(f"MODEL{'2'.rjust(9, ' ')}\n")
        atom_count = 1
        subchain = fitted_protein_2_chain
        for r in fitted_protein_2:
            res_name = r.name.rjust(3, " ")
            res_num = str(r.seqid.num).rjust(4, " ")
            for a in r:
                atom_num = str(atom_count).rjust(5, " ")
                atom_name = a.name.ljust(4, " ")
                x = str(round(a.pos.x, 3)).rjust(8, " ")
                y = str(round(a.pos.y, 3)).rjust(8, " ")
                z = str(round(a.pos.z, 3)).rjust(8, " ")
                row = f"ATOM  {atom_num} {atom_name} {res_name} {subchain}{res_num}    {x}{y}{z}\n"
                fw.write(row)
                atom_count += 1
        fw.write("ENDMDL\n")
    except Exception as e:
        traceback.print_exc()
        print(e)


def write_final_output_pml(folder_name: str, protein, domains, fixed_domain_id, bending_residues):
    bend_res_colour = "[0  ,255,0  ]"
    dom_colours = ["[0  ,0  ,255]", "[255,0  ,0  ]", "[255,255,0  ]", "[255,100,255]", "[0  ,255,255]"]
    print(domains)
    for b in bending_residues.values():
        for k, g in groupby(enumerate(b), lambda ix: ix[0] - ix[1]):
            print(list(map(itemgetter(1), g)))

    try:
        dir_path_str = f"{output_path}/{folder_name}"
        dir_path = Path(dir_path_str)
        if not dir_path.exists():
            dir_path.mkdir(parents=True)
        fw = open(f"{dir_path}/{folder_name}.pml", "w")
        fw.write("reinitialize\n")
        fw.write(f"load {folder_name}.pdb\n")
        fw.write(f"bg_color white\n")
        fw.write("color grey\n")

        fixed_dom_segments = domains[fixed_domain_id].segments
        util_res = protein.utilised_residues_indices

        fixed_dom_res_reg = []
        all_bend_res_indices = []
        for b in bending_residues.values():
            all_bend_res_indices.extend(b)
        for s in range(fixed_dom_segments.shape[0]):
            reg = []
            for i in range(fixed_dom_segments[s][0], fixed_dom_segments[s][1]+1):
                index = util_res[i]
                if index not in all_bend_res_indices:
                    reg.append(index)
            fixed_dom_res_reg.extend(group_continuous_regions(reg))

        for s in range(len(fixed_dom_res_reg)):
            if s == 0:
                sel_reg_str = f"select region0, resi {fixed_dom_res_reg[s][0]}-{fixed_dom_segments[s][1]}\n"
            else:
                sel_reg_str = f"select region0, region0 + resi {fixed_dom_segments[s][0]}-{fixed_dom_segments[s][1]}\n"
            fw.write(sel_reg_str)
        fw.write(f"set_color colour0 = {dom_colours[0]}\n")
        fw.write("color colour0, region0\n")

        region_count = 1
        for domain in domains:
            dyn_dom_res_reg = []
            if domain.domain_id == fixed_domain_id:
                continue
            segments = domain.segments
            dom_bend_res = bending_residues[domain.domain_id]
            for s in range(segments.shape[0]):
                reg = []
                for i in range(segments[s][0], segments[s][1]+1):
                    index = util_res[i]
                    if index not in dom_bend_res:
                        reg.append(index)

                dyn_dom_res_reg.extend(group_continuous_regions(reg))

            for s in range(len(dyn_dom_res_reg)):
                if s == 0:
                    sel_reg_str = f"select region{region_count}, resi {dyn_dom_res_reg[s][0]}-{dyn_dom_res_reg[s][1]}\n"
                else:
                    sel_reg_str = f"select region{region_count}, region{region_count} + resi {dyn_dom_res_reg[s][0]}-{dyn_dom_res_reg[s][1]}\n"
                fw.write(sel_reg_str)
            fw.write(f"set_color colour{region_count} = {dom_colours[region_count]}\n")
            fw.write(f"color colour{region_count}, region{region_count}\n")

            region_count += 1

        bend_res_groups = group_continuous_regions(all_bend_res_indices)
        for g in bend_res_groups:
            fw.write(f"select region{region_count}, resi {g[0]}-{g[1]}\n")
            fw.write(f"set_color colour{region_count} = {bend_res_colour}\n")
            fw.write(f"color colour{region_count}, region{region_count}\n")
            region_count += 1

        fw.write("set dash_gap, 0\n")
        fw.write("set dash_radius, 0.2\n")

    except Exception as e:
        traceback.print_exc()
        print(e)


def write_w5_info_file(protein_1_name: str, protein_2_name: str, param: dict, domains: list,
                       fixed_domain_id: int):
    try:
        fw = open(f"output/w5_info/{protein_1_name}{param['chain1id']}_{protein_2_name}{param['chain2id']}.w5_info", "w")
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
        for domain in domains:
            if domain.domain_id == fixed_domain_id:
                fw.write("FIXED DOMAIN\n")
                fw.write(f"DOMAIN NUMBER: \t {fixed_domain_id} (coloured yellow for rasmol)\n")
                residue_str = ""
                for s in range(domain.segments.shape[0]):
                    if s == domain.segments.shape[0] - 1:
                        residue_str += str(domain.segments[s][0]) + " - " + str(domain.segments[s][1])
                    else:
                        residue_str += str(domain.segments[s][0]) + " - " + str(domain.segments[s][1]) + " , "
                fw.write(f"RESIDUE NUMBERS: \t{residue_str}\n")
                fw.write(f"SIZE: \t{domain.num_residues}\n")
                fw.write(f"BACKBONE RMSD ON THIS DOMAIN: \t{round(domain.rmsd, 3)}A\n")
                break

        domain_count = 1
        for domain in domains:
            if domain.domain_id != fixed_domain_id:
                fw.write("------------------------------------------------------------------------------\n")
                fw.write(f"MOVING DOMAIN (RELATIVE TO FIXED DOMAIN),  PAIR {domain_count}\n")
                domain_count += 1
                colour_str = ""
                if domain_count == 1:
                    colour_str = "blue"
                elif domain_count == 2:
                    colour_str = "green"
                elif domain_count == 3:
                    colour_str = "red"
                else:
                    colour_str = "purple"
                fw.write(f"DOMAIN NUMBER: \t {domain.domain_id} (coloured {colour_str} for rasmol)\n")
                residue_str = ""
                for s in range(domain.segments.shape[0]):
                    if s == domain.segments.shape[0] - 1:
                        residue_str += str(domain.segments[s][0]) + " - " + str(domain.segments[s][1])
                    else:
                        residue_str += str(domain.segments[s][0]) + " - " + str(domain.segments[s][1]) + " , "
                fw.write(f"RESIDUE NUMBERS: \t{residue_str}\n")
                fw.write(f"SIZE: \t{domain.num_residues}\n")
                fw.write(f"BACKBONE RMSD ON THIS DOMAIN: \t{round(domain.rmsd, 3)}A\n")

    except Exception as e:
        print(e)
        return False
    return True


def group_continuous_regions(data: list):
    groups = []
    for k, g in groupby(enumerate(data), lambda ix: ix[0] - ix[1]):
        temp = list(map(itemgetter(1), g))
        groups.append([temp[0], temp[-1]+1])
    return groups


