import math
from receptor import Receptor
from results import Results
from itertools import permutations

def find_attributes(receptor, results):
    """
    Identify pharmacophores from the docking results and returns the top score pharmacophores
    :param receptor: the receptor in .pdbqt file
    :param results: the result of AutoDock Vina in .pbdqt file
    :return: top score pharmacophores
    """
    rec_atoms = receptor.atoms
    lig_atoms = results.atoms
    attribute = []
    unique_pairs = []

    for atom in rec_atoms:
        atom_type = atom.autodock_atom_type
        res_name = atom.residue_name

        # 1. Identify ion - dipole
        if res_name in ("ARG", "LYS") and atom.autodock_atom_type in ("N"):
            # positive charge from the rec and partial negative charge from the ligand
            for l_atom in lig_atoms:
                if l_atom.atom_name not in ("C"):
                    if l_atom.partial_charge < -0.010:
                        distance = math.sqrt(
                            (l_atom.x_coord - atom.x_coord) ** 2 + (l_atom.y_coord - atom.y_coord) ** 2 +
                            (l_atom.z_coord - atom.z_coord) ** 2)
                        if 1.4< distance < 5.6:
                            unique_pairs.append((atom, l_atom, distance))
                            attribute.append(("pos ion - neg par", distance, atom, l_atom))

        if res_name in ("ASP", " GLU"):
            # TODO check the atom partial charge threshold
            if l_atom.atom_name not in ("C") and atom.autodock_atom_type in ("OA"):
                if l_atom.partial_charge > +0.010 and atom.partial_charge < -0.3:
                    distance = math.sqrt(
                        (l_atom.x_coord - atom.x_coord) ** 2 + (l_atom.y_coord - atom.y_coord) ** 2 +
                        (l_atom.z_coord - atom.z_coord) ** 2)
                    if 1.4< distance < 5.6:
                        unique_pairs.append((atom, l_atom, distance))
                        attribute.append(("neg ion - pos par", distance, atom, l_atom))

        # 2. Identify hydrogen bonding
        # if the atom is hydrogen bond acceptor, find hydrogen bond donor
        if atom_type in ("NA", "OA", "F", "S"):
            for l_atom in lig_atoms:
                # print(" LIGAND: " + l_atom.autodock_atom_type)
                if l_atom.autodock_atom_type in ("HD"):
                    # Calculate 3D distance
                    distance = math.sqrt((l_atom.x_coord - atom.x_coord)**2 + (l_atom.y_coord - atom.y_coord)**2 +
                                         (l_atom.z_coord - atom.z_coord)**2)
                    if 2.1< distance < 3.9:
                        new_atribute = (atom, l_atom, distance)
                        if new_atribute not in unique_pairs:
                            attribute.append(("hydrogen bonding", distance, atom, l_atom))
        # if the atom is hydrogen bond donor, find hydrogen bond acceptor
        if atom_type in ("HD"):
            for l_atom in lig_atoms:
                if l_atom.autodock_atom_type in ("NA", "OA", "F", "S"):
                    # Calculate 3D distance
                    distance = math.sqrt((l_atom.x_coord - atom.x_coord) ** 2 + (l_atom.y_coord - atom.y_coord) ** 2 +
                                         (l_atom.z_coord - atom.z_coord) ** 2)
                    if 2.1< distance < 3.9:
                        new_atribute = (atom, l_atom, distance)
                        if new_atribute not in unique_pairs:
                            attribute.append(("hydrogen bonding", distance, atom, l_atom))

        # 3. Identify dipole - dipole interaction
        # TODO should i include all the atom in the receptor (also the backbone interaction of other amino acid)?
        if res_name in ("SER", "THR", "CYS", "ASN", "GLN") and atom.autodock_atom_type in ("OA", "S", "N"):
            # positive charge from the rec and partial negative charge from the ligand
            for l_atom in lig_atoms:
                if l_atom.partial_charge > +0.010 and atom.partial_charge < -0.010:
                    distance = math.sqrt(
                            (l_atom.x_coord - atom.x_coord) ** 2 + (l_atom.y_coord - atom.y_coord) ** 2 +
                            (l_atom.z_coord - atom.z_coord) ** 2)
                    if 1.4 < distance < 5.6:
                        new_atribute = (atom, l_atom, distance)
                        if new_atribute not in unique_pairs:
                            attribute.append(("neg par - pos par", distance, atom, l_atom))

        # 4. Identify hydrophobic interaction
        if res_name in ("ALA", "VAL", "ILE", "LEU", "PHE", "MET", "TYR", "TRP"):
            for l_atom in lig_atoms:
                if l_atom.autodock_atom_type in ("A ") and atom.autodock_atom_type in ("C ", "C"):
                    # Calculate 3D distance
                    distance = math.sqrt((l_atom.x_coord - atom.x_coord) ** 2 + (l_atom.y_coord - atom.y_coord) ** 2 +
                                         (l_atom.z_coord - atom.z_coord) ** 2)
                    if 1.0 < distance < 5.9:
                        new_atribute = (atom, l_atom, distance)
                        if new_atribute not in unique_pairs:
                            attribute.append(("ldf", distance, atom, l_atom))
    return attribute

def get_unique_attributes(attribute):
    """

    :param attribute: a list of all the attribute
    :return: list of unique interaction from all the poses with the top distance
    """
    unique_pairs = {}
    for a in attribute:
        new_pair = a[2].residue_name + a[2].residue_sequence, a[2].atom_name,  a[3].atom_name, a[3].serial_number, a[0]
        # print(new_pair)
        if new_pair not in unique_pairs.keys():
            unique_pairs.update({new_pair : a[1]})
        else:
            old_distance = unique_pairs.get(new_pair)
            if old_distance < a[1]:
                unique_pairs.update({new_pair: a[1]})
    return unique_pairs


def generate_pharmacophore (unique_pairs, num_attrib):
    """
    generate pharmacophores based on number of attribute
    """
    pharmacophores = []
    # Generate different combinations of pharmacophore from (unique_pairs)
    # and each combination has num_attrib items.
    perm = permutations(unique_pairs, num_attrib)
    unique_comb = []
    for i in list(perm):
        residue_name = []
        for j in i:
            new_residue = (j[0])
            # print(new_residue)
            if new_residue not in residue_name:
                residue_name.append(j[0])
        if len(residue_name) == 3:
            score = scoring_function(i)
            interaction = {}
            interaction.update({"bond": i})
            interaction.update({"score": score})
            pharmacophores.append(interaction)

    pharmacophores.sort(key=getScore, reverse=True)
    return pharmacophores

def getScore (pharmacophore):
    """
    A function that return the score value
    """
    return pharmacophore["score"]

def scoring_function(list):
    """
    Return the score for the pharmacophore
    """
    score = 0;
    for interaction in list:
        # distance-dependent function
        bonding_interaction = interaction[4]
        if bonding_interaction in ("pos ion - neg par", "neg ion - pos par"):
           score += 4
        elif bonding_interaction in ("hydrogen bonding"):
            score += 3
        elif bonding_interaction in ("neg par - pos par"):
            score += 2
        elif bonding_interaction in ("ldf"):
            score += 1
    return score


if __name__ == '__main__':
    center = [51.512, 22.11, 104.839]
    size = [30, 30, 30]
    rec_path = "/Users/mythanhthaonguyen/PycharmProjects/Thesis_DockingAnalysis/data/7jtl_mol.pdbqt"
    receptor = Receptor(rec_path, center, size)
    results_path = "/Users/mythanhthaonguyen/PycharmProjects/Thesis_DockingAnalysis/data/ZINC857.pdbqt"
    docking_result = Results(results_path, center, size)
    attributes = find_attributes(receptor, docking_result)
    # for a in attributes:
        # print(a[0] + " " + str(a[1]) + " between: " + a[2].residue_name + " and " + a[3].atom_name + " " + a[3].serial_number)
        # print(a)
    dict = get_unique_attributes(attributes)
    # # for pair in dict:
    # #     print(pair)
    # #     print(dict.get(pair))
    top_pharma = generate_pharmacophore(dict, 3)
    for i in range(10):
        print(top_pharma[i])


