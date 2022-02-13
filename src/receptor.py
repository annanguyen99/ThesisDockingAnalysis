import pathlib
from atom import Atom

class Receptor:
    """
    Class used for reading receptor from .pdbqt file type and save all the atoms that in the binding pocket
    """
    def __init__(self, path, center, size):
        # All the atoms that within the docking box
        self.atoms = []
        # All the unique residues that within the docking box
        self.residues_id = []
        # The center of the box
        self.center = center
        # The size of the box
        self.size = size

        # Read the receptor file
        with open(path) as f:
            # Read in the line and create new atoms
            # until the end of the file
            for line in f:
                if line[:6] in ("ATOM  ", "HETATM"):
                    new_atom = Atom(line)
                    # Check to see if the atom is within the box
                    # if (float(new_atom.x_coord) < (center[0] + size[0]) and float(new_atom.x_coord) > (center[0] - size[0]))and \
                    #         (float(new_atom.y_coord) < (center[1] + size[1]) and float(new_atom.x_coord) > (center[1] - size[1])) and \
                    #         (float(new_atom.z_coord) < (center[2] + size[2]) and float(new_atom.z_coord) > (center[2] - size[2])):

                    # If the atom is not water, add to the list
                    if str(new_atom.residue_name) != "HOH":
                        self.atoms.append(new_atom)
                        # Add the unique residue and residue id to the list (and not water)
                        if (new_atom.residue_name, new_atom.residue_sequence) not in self.residues_id:
                            self.residues_id.append((new_atom.residue_name, new_atom.residue_sequence))

if __name__ == "__main__":

    # Quick testing

    # Project path
    project_root = pathlib.Path(__file__).parent.parent
    center = [51.512, 22.11, 104.839]
    size = [30, 30, 30]
    receptor_path = project_root / 'data' / '7jtl.pdbqt'
    receptor = Receptor(receptor_path, center, size)

    # Print all the atoms' name and its x-coordinate
    for atom in receptor.atoms:
        print(atom.atom_name + ": " + str(atom.x_coord))
