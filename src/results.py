import pathlib
from atom import Atom

class Results:
    """
    Class used for reading AutoDock Vina docking results and save the atoms of the ligands
    """
    def __init__(self, path, center, size):
        self.num_pose = 0
        self.atoms = []
        with open(path) as f:
            for line in f:
                if line[:6] in ("MODEL "):
                    self.num_pose += 1
                if line[:6] in ("ATOM  ", "HETATM", "ATOM"):
                    new_atom = Atom(line)
                    # if (float(new_atom.x_coord) < (center[0] + size[0]) and float(new_atom.x_coord) > (center[0] - size[0]))and \
                    #         (float(new_atom.y_coord) < (center[1] + size[1]) and float(new_atom.x_coord) > (center[1] - size[1])) and \
                    #         (float(new_atom.z_coord) < (center[2] + size[2]) and float(new_atom.z_coord) > (center[2] - size[2])):
                    self.atoms.append(new_atom)

if __name__ == "__main__":

    # Quick testing
    center = [51.512, 22.11, 104.839]
    size = [34, 34, 34]

    project_root = pathlib.Path(__file__).parent.parent
    results_path = project_root / 'data' / 'ZINC596.pdbqt'
    docking_result = Results(results_path, center, size)

    # Print out all the atoms
    for atom in docking_result.atoms:
        print(atom.atom_name + ": " + str(atom.serial_number))

