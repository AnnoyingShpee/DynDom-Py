import sys
import gemmi


# file = "pdb_files/1lfg.cif"
file = "pdb_files/1lfg.pdb"

structure: gemmi.Structure = gemmi.read_structure(file)
unit_cell: gemmi.Structure.cell = structure.cell
space_group: gemmi.Structure.spacegroup_hm = structure.spacegroup_hm
print(f"Structure =", structure)
print(f"Unit Cell =", unit_cell)
print(f"Space Group =", space_group)

# Iterate over atoms in the structure
for model in structure:  # gemmi.Model
    # print(model)  # Prints <gemmi.Model 1 with 3 chain(s)>
    for chain in model:  # gemmi.Chain
        # print(chain)  # Prints 3 <gemmi.Chain * with ** res>
        for residue in chain:  # gemmi.Residue
            # print(residue)
            # print(type(residue))
            for atom in residue:  # gemmi.Atom
                # atomic: gemmi.Atom = atom
                # print(atomic.pos)
                print(atom)
