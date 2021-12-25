class Residue:
    """
    Create an residue object (might not need this class)
    """
    def __init__(self, name, id):
        # The abbreviation of the residue
        self.name = name
        self.id = id
        self.hba = False
        self.hbd = False
        self.pos_ion = False
        self.neg_ion = False
        self.hydrophobic = False
        self.list = ["ALA, CYS, ASP, GLU, PHE, GLY, HIS, ILE, LYS, LEU, MET, ASN,"
                     "PRO, GLN, ARG, SER, THR, VAL, TRP, TYR"]

        # TODO check to see if histidine is positive at pH 7.4
        if self.name in ("ARG", "LYS"):
            self.pos_ion = True

        if self.name in ("ASP", " GLU"):
            self.neg_ion = True

        if self.name in ("ALA", "VAL", "ILE", "LEU", "PHE", "MET", "TYR", "TRP"):
            self.hydrophobic = True

        if self.name in ("ARG", "HIS", "ASP", "GLU", "SER", "THR", "ASN", "GLN", "CYS", "MET", "TYR"):
            self.hba = True

        if self.name in ("ARG", "HIS", "LYS", "SER", "THR", "ASN", "GLN", "CYS", "TYR"):
            self.hbd = True



