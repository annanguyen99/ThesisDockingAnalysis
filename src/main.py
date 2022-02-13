import math
import pathlib
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

    # Go through all the atoms in the receptor
    for atom in rec_atoms:

        # Save the atom type
        atom_type = atom.autodock_atom_type

        # Save the residue of the atom
        res_name = atom.residue_name

        # 1. Identify ion - dipole
        # positive charge from the receptor and partial negative charge from the ligand
        if res_name in ("ARG", "LYS") and atom.autodock_atom_type in ("N"):
            for l_atom in lig_atoms:
                if l_atom.atom_name not in ("C"):
                    if l_atom.partial_charge < -0.010:
                        # Calculate 3D distance
                        distance = math.sqrt(
                            (l_atom.x_coord - atom.x_coord) ** 2 + (l_atom.y_coord - atom.y_coord) ** 2 +
                            (l_atom.z_coord - atom.z_coord) ** 2)
                        if 1.4< distance < 5.6:
                            unique_pairs.append((atom, l_atom, distance))
                            attribute.append(("pos ion - neg par", distance, atom, l_atom))

        # negative charge from the receptor and partial positive charge from the ligand
        if res_name in ("ASP", " GLU"):
            if l_atom.atom_name not in ("C") and atom.autodock_atom_type in ("OA"):
                if l_atom.partial_charge > +0.010 and atom.partial_charge < -0.3:
                    # Calculate 3D distance
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
    # Return a list of all receptor's atom and ligand's atom interaction
    return attribute

def get_unique_attributes(attribute):
    """
    This method gets all the unique interactions between the receptor and ligand
    :param attribute: a list of all the attribute
    :return: list of unique interaction from all the poses with the top distance
    """
    # Contain all the unique interaction between receptor's atoms and ligand's atoms
    unique_pairs = {}
    # Contain all the unique residue that have interaction from receptor
    unique_residues = {}

    # Go through all the attributes
    for a in attribute:
        # Residue id includes: the receptor's residue and its unique residue sequence number
        # eg. LYS47
        residue_id = a[2].residue_name + a[2].residue_sequence
        # New pair includes: the receptor's residue, the receptor's atom name,
        # the ligand's atom name, and its serial number, and the interaction's name
        new_pair = residue_id, a[2].atom_name,  a[3].atom_name, a[3].serial_number, a[0]

        # Check to see if the interaction is already in the unique_pairs
        if new_pair in unique_pairs.keys():
            # Calculate the new distance if it's already in the unique interactions map
            old_distance = unique_pairs.get(new_pair)
            if old_distance < a[1]:
                # If the new distance is better than the old distance, replace the distance
                # in the map with the new distance
                unique_pairs.update({new_pair: a[1]})

        # If not the exactly rec's atom -  lig's atom interaction in the map,
        # check to see if the rec's residue already in the unique residue map
        elif (residue_id) in unique_residues.keys():
            # If it is in the map, compare the distance of the old distance and the
            # new distance
            old_pairs = unique_residues.get(residue_id)
            old_distance = unique_pairs.get(old_pairs)
            if old_distance < a[1]:
                # If the new distance is better than the old distance, replace the distance
                # in the map with the new distance
                unique_pairs.pop(old_pairs)
                unique_pairs.update({new_pair: a[1]})
                unique_residues.update({(residue_id) : new_pair})
        else:
        # If the interaction is in neither list, add to both of the list
            unique_pairs.update({new_pair : a[1]})
            unique_residues.update({(residue_id) : new_pair})

    return unique_pairs


def generate_pharmacophore (unique_pairs, num_attrib):
    """
    This method generates the generate different combination of the pharmacophores from
    the unique interactions of receptor and ligands
    and sorted the list
    :param unique_pairs: the unique interactions in the result file
    :param num_attrib: the ideal number of attributes in the pharmacophore
    :return: the list of pharmacophores
    """
    pharmacophores = []
    perm = permutations(unique_pairs, num_attrib)

    for i in perm:
        residue_name = []
        for j in i:
            new_residue = (j[0])
            if new_residue not in residue_name:
                residue_name.append(j[0])
        if len(residue_name) == num_attrib:
            score = scoring_function(i)
            interaction = {}
            interaction.update({"bond": i})
            interaction.update({"score": score})
            pharmacophores.append(interaction)

    pharmacophores.sort(key=getScore, reverse=True)
    return pharmacophores

def getScore (pharmacophore):
    """
    A get function that return the score value
    """
    return pharmacophore["score"]

def scoring_function(list):
    """
    Return the score for the pharmacophore
    :param list: the list of attributes in a single pharmacophore
    :return: the total score of one pharmacophore based on its intermolecular force s
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
    project_root = pathlib.Path(__file__).parent.parent
    receptor_path = project_root / 'data' / '7jtl.pdbqt'
    results_path = project_root / 'data' / 'ZINC596.pdbqt'
    print()
    print("The receptor path: " + str(receptor_path))
    print("The result path: " + str(results_path))

    # Center of the docking box
    center = [51.512, 22.11, 104.839]
    # The size of the docking box
    size = [38, 48, 34]
    # Create the receptor object
    receptor = Receptor(receptor_path, center, size)
    # Create an results object based on the result file
    docking_result = Results(results_path, center, size)
    attributes = find_attributes(receptor, docking_result)

    # Get all the unique interactions
    dict = get_unique_attributes(attributes)

    # Generate the pharmacophores
    num_attribute = 3
    top_pharma = generate_pharmacophore(dict, num_attribute)
    print()
    print("This is a list of top 10 scored potential pharmacophores with " + str(num_attribute) + " attribute per pharmacophore")
    # Print out top 10 scored pharmacophores
    for i in range(10):
        pharma = top_pharma[i]
        print()
        print("Pharmacophore " + str(i + 1) + " with score = " + str(pharma.get("score")) )
        interactions = pharma.get("bond")
        for j in range(len(interactions)):
            interaction = interactions[j]
            print("Receptor's atom " + interaction[1] + " of residue " + interaction[0] +
                  " and ligand's atom " + interaction[2] +interaction[3])
