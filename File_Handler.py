

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

