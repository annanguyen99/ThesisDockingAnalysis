
class Atom:
    """
    Class used for reading one line of ATOM in PDBQT file type and create an ATOM object
    """

    def __init__(self, line: str):
        self.record_name = line[:6]
        self.serial_number = line[9:11]
        self.atom_name = line[13:16]
        self.residue_name = line[17:20]
        self.chainID = line[21]
        self.residue_sequence = line[24:26]
        self.x_coord = float(line[32:38])
        self.y_coord = float(line[40:46])
        self.z_coord = float(line[48:54])
        self.partial_charge = float(line[70:76])
        self.autodock_atom_type = line[77:79]


if __name__ == "__main__":

    # Quick testing
    # atom = Atom("ATOM      1  N   GLN A  18      40.917  33.173 126.601  1.00 73.92    -0.058 N ")
    atom = Atom("ATOM      6  C   UNL     1      48.774  29.090  98.883  0.00  0.00    +0.000 A ")
    print(atom.record_name)
    print(atom.serial_number)
    print(atom.atom_name)
    print(atom.residue_name)
    print(atom.chainID)
    print(atom.residue_sequence)
    print(atom.x_coord)
    print(atom.y_coord)
    print(atom.z_coord)
    print(atom.partial_charge)
    print(atom.autodock_atom_type)