import sys
import gemmi

class Protein_Structure:
    def __init__(self, file_name: str):
        self.file_name: str = file_name
        self.structure: gemmi.Structure = ...
        self.structure_metadata = None
        self.file_type: str = ""  # Can be pdb or cif file
        self.check_file_format(file_name)

    def check_file_format(self, file_name: str) -> None:
        if file_name.endswith(".pdb"):
            self.file_type = "pdb"
        elif file_name.endswith(".cif"):
            self.file_type = "cif"
        else:
            print("Unsupported file type")
            return

    # greeted = set()
    # for path in sys.argv[1:]:
    #     try:
    #         doc = gemmi.cif.read_file(path)  # copy all the data from cif file
    #         block = doc.sole_block()  # cif has exactly one block
    #         for element in block.find_loop("_atom_site.type_symbol"):
    #             if element not in greeted:
    #                 print("Hello " + element)
    #                 greeted.add(element)
    #     except Exception as e:
    #         print("Oops. %s" % e)
    #         sys.exit(1)

    def get_cif_document(self):
        self.structure: gemmi.cif.Document = gemmi.cif.read_file(self.file_name)
        print(self.structure.as_string())
        # self.structure = gemmi.read_structure(self.file_name, format=gemmi.CoorFormat.Mmcif)
        # Returns a gemmi.cif.Block object. Can
        # self.structure_metadata: gemmi.cif.Block = self.structure.make_mmcif_headers()

    def get_pdb_structure(self):
        self.structure = gemmi.read_structure(self.file_name, format=gemmi.CoorFormat.Pdb)
        # Returns a string of the whole structure's details
        self.structure_metadata = self.structure.make_pdb_headers()

    def get_structure(self):
        self.structure = gemmi.read_structure(self.file_name, format=gemmi.CoorFormat.Detect)
        self.structure_metadata = self.structure

    def print_pipeline(self):
        if self.file_type == "pdb":
            self.get_pdb_structure()
        elif self.file_type == "cif":
            self.get_cif_document()

        # print(f"Structure type = {type(self.structure)}")
        # print(self.structure)
        # print(f"Detailed type = {type(self.structure_metadata)}")
        # print(self.structure_metadata)

    # https://gemmi.readthedocs.io/en/latest/mol.html#mcra
    # def iterate_hierarchy_levels(self):


file = "data/cif/1lfg.cif"
# file = "pdb/1lfg.pdb"
test = Protein_Structure(file)
test.print_pipeline()

