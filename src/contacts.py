"""
Author: Rafael Lemos - rafaellemos42@gmail.com
Date: 12/08/2024

License: MIT License
"""

from math import dist
from numpy import dot, arccos, degrees
from numpy.linalg import norm
from copy import deepcopy
from collections import defaultdict

from src.classes import Contact
from src.distances import distances
import src.conditions as conditions


def contact_detection(protein, region, chains, interface, custom_distances, epsilon, uncertainty_flags, local_contact_types):
    """
    Detects contacts between atoms in a given protein.

    Args:
        protein (Protein): The protein object containing chain, residue and atom objects.

    Returns:
        list: A list of Contact objects representing the detected contacts.
    """

    residues = list(protein.get_residues())
    contacts = []
    interface_res = []
    uncertain_contacts = []
    max_ca_distance = 20.47 # 0.01 higher than the Arg-Arg pair

    contacts_by_pair = defaultdict(list) # empty list under the key if it doesn't exist
    contact_registry = defaultdict(set) # empty set under the key if it doesn't exist

    categories = custom_distances if custom_distances else conditions.categories
    if epsilon > 0:
        max_ca_distance += epsilon
        updated_distances = {key: value + epsilon for key, value in distances.items()}
    else:
        updated_distances = distances
    
    if interface:
        with open(interface,"r") as f:
            for line in f:
                interface_res.append(line.strip())
        
    for i, residue1 in enumerate(residues[1:]):
        for _, residue2 in enumerate(residues[i+1:], start=i+1):
            
            if residue1.resnum == residue2.resnum and residue1.chain.id == residue2.chain.id: # ignores same residue
                continue
            
            if region and (residue1.resnum not in region or residue2.resnum not in region):
                continue

            if chains and (residue1.chain.id not in chains or residue2.chain.id not in chains):
                continue

            if len(residue1.atoms) > 1 and len(residue2.atoms) > 1:
                ca1, ca2 = residue1.atoms[1], residue2.atoms[1] # alpha carbons

                distance_ca = dist((ca1.x, ca1.y, ca1.z), (ca2.x, ca2.y, ca2.z))
                
                # filter distant residues (static value then specific values)
                if distance_ca > max_ca_distance:
                    continue
                else:
                    key = ''.join(sorted((residue1.resname, residue2.resname)))
                    if distance_ca > (updated_distances[key] + epsilon):
                        continue

            else:
                continue

            pair_key = (residue1.chain.id, residue1.resnum, residue2.chain.id, residue2.resnum)

            # CHECKING FOR AROMATIC STACKINGS
            if residue1.ring and residue2.ring:
                ring1, ring2 = residue1.atoms[-1], residue2.atoms[-1] # RNG atoms
                if interface and ring1.entity == ring2.entity:
                    continue
                
                distance = dist((ring1.x, ring1.y, ring1.z), (ring2.x, ring2.y, ring2.z))
                angle = calc_angle(residue1.normal_vector, residue2.normal_vector)
                
                aromatic_range = categories['aromatic']
                if aromatic_range[0] <= distance <= aromatic_range[1]:
                    if (160 <= angle < 180) or (0 <= angle < 20):
                        stack_type = "-parallel"
                    elif (80 <= angle < 100):
                        stack_type = "-perpendicular"
                    else:
                        stack_type = "-other"

                    contacts_by_pair[pair_key].append({
                        'protein_id': protein.id,
                        'chain1': residue1.chain.id,
                        'resnum1': residue1.resnum,
                        'resname1': residue1.resname,
                        'atomname1': ring1.atomname,
                        'chain2': residue2.chain.id,
                        'resnum2': residue2.resnum,
                        'resname2': residue2.resname,
                        'atomname2': ring2.atomname,
                        'distance': float(f"{distance:.2f}"),
                        'type': "stacking"+stack_type,
                        'atom1': ring1,
                        'atom2': ring2,
                    })
                    
            for atom1 in residue1.atoms:
                for atom2 in residue2.atoms:
                    
                    if interface:
                        residue_interface_key = f"{residue1.chain.id},{residue1.resnum},{residue1.resname}"
                        if (atom1.entity == atom2.entity) or (residue_interface_key not in interface_res):
                            continue
                    
                    name1 = f"{atom1.residue.resname}:{atom1.atomname}" # matches the pattern from conditions dictionary
                    name2 = f"{atom2.residue.resname}:{atom2.atomname}"

                    if name1 in local_contact_types and name2 in local_contact_types: # excludes the RNG atom and any different other

                        distance = dist((atom1.x, atom1.y, atom1.z), (atom2.x, atom2.y, atom2.z))
                        
                        if distance <= 6: # max distance for contacts

                            for contact_type, distance_range in categories.items():

                                if contact_type == 'hydrogen_bond' and (abs(residue2.resnum - residue1.resnum) <= 3): # skips alpha-helix for h-bonds
                                    continue

                                if distance_range[0] <= distance <= distance_range[1]: # fits the range

                                    def get_props(name, contact_type = contact_type):
                                        if name in uncertainty_flags and contact_type in ['attractive','repulsive','salt_bridge']:
                                            return resolve_uncertainty(name, uncertainty_flags, local_contact_types)
                                        if contact_type == 'disulfide_bond':
                                            return name
                                        return local_contact_types[name]

                                    props1 = get_props(name1)
                                    props2 = get_props(name2)

                                    if not conditions.contact_conditions[contact_type](props1, props2):
                                        continue
                                    
                                    stored_types = contact_registry[pair_key] # empty set if pair_key is new
                                                                                            
                                    if contact_type == 'salt_bridge':
                                        if 'attractive' in stored_types:
                                            # filters attractives out in-place (salt_bridge has priority)
                                            contacts_by_pair[pair_key] = [c for c in contacts_by_pair[pair_key] if c['type'] != 'attractive']
                                            stored_types.remove('attractive')
                                        if 'salt_bridge' in stored_types:
                                            continue # skip adding duplicate
                                        stored_types.add('salt_bridge')
                                    elif contact_type == 'attractive':
                                        if 'attractive' in stored_types or 'salt_bridge' in stored_types:
                                            continue
                                        stored_types.add('attractive')
                                    elif contact_type == 'repulsive':
                                        if 'repulsive' in stored_types:
                                            continue
                                        stored_types.add('repulsive')
                                    else:
                                        stored_types.add(contact_type)
                                    
                                    contacts_by_pair[pair_key].append({
                                        'protein_id': protein.id,
                                        'chain1': residue1.chain.id,
                                        'resnum1': residue1.resnum,
                                        'resname1': residue1.resname,
                                        'atomname1': atom1.atomname,
                                        'chain2': residue2.chain.id,
                                        'resnum2': residue2.resnum,
                                        'resname2': residue2.resname,
                                        'atomname2': atom2.atomname,
                                        'distance': float(f"{distance:.2f}"),
                                        'type': contact_type,
                                        'atom1': atom1,
                                        'atom2': atom2,
                                        'is_uncertain':(name1 in uncertainty_flags or name2 in uncertainty_flags) and contact_type in ['attractive', 'repulsive', 'salt_bridge'],
                                    })
    
    contacts, count_types = create_contacts(contacts_by_pair)
                                            
    return contacts, interface_res, count_types, uncertain_contacts


def show_contacts(contacts):
    """
    Formats and summarizes contact information to be outputted to a file. Only works with the -o flag.

    Args:
        contacts (list): A list of Contact objects of a given protein.

    Returns:
        str: A formatted string summarizing the contact information.
    """
    
    output = []
    
    output.append("Chain1,Res1,ResName1,Atom1,Chain2,Res2,ResName2,Atom2,Distance,Type")
    for contact in contacts:
        output.append(contact.print_text())
        
    return "\n".join(output) # returns as a string to be written directly into the file


def calc_angle(vector1, vector2):
    """
    Calculates the angle between two ring vectors of aromatic residues

    Args:
        vector1 (tuple): The first vector (x, y, z).
        vector2 (tuple): The second vector (x, y, z).

    Returns:
        float: The angle between the vectors in degrees.
    """
    
    dot_product = dot(vector1, vector2)
    magnitude_product = norm(vector1) * norm(vector2) # normalizes the dot product
    angle = arccos(dot_product / magnitude_product) # angle in radians   
    
    return degrees(angle)


def change_protonation(ph, silent):
    from src.process import log
    
    pka_table = {
        'R': 12.48,
        'K': 10.79,
        'H': 6.04,
        'D': 3.86,
        'E': 4.25,
        'C': 8.33,
        'Y': 10.07,
    }
    
    pH_sensitive_atoms = {
        'R': ['NE', 'CZ', 'NH1', 'NH2'],
        'K': ['NZ'],
        'H': ['ND1', 'NE2'],
        'D': ['OD1', 'OD2'],
        'E': ['OE1', 'OE2'],
        'C': ['SG'],
        'Y': ['OH'],
    }
    
    uncertainty_flags = {}
    local_contact_types = deepcopy(conditions.contact_types)
    
    for key, value in local_contact_types.items():
        resname, atomname = key.split(":")
        if resname in pka_table and atomname in pH_sensitive_atoms.get(resname, []):
            pka = pka_table[resname]
            delta = abs(ph - pka)
            
            original_pos = value[2]
            original_neg = value[3]
            
            new_pos, new_neg = original_pos, original_neg  # Default: no change            
            
            if resname in ['D', 'E', 'C', 'Y']:  # Acidic
                if delta < 2.0:
                    new_pos = 0
                    new_neg = 0
                    uncertainty_flags[key] = {'neg': True}
                else:
                    is_deprotonated = ph > pka
                    new_pos = 0
                    new_neg = 1 if is_deprotonated else 0

            elif resname in ['R', 'K', 'H']:  # Basic
                if delta < 2.0:
                    new_pos = 0
                    new_neg = 0
                    uncertainty_flags[key] = {'pos': True}
                else:
                    is_protonated = ph < pka
                    new_pos = 1 if is_protonated else 0
                    new_neg = 0
                        
            if (original_pos != new_pos) or (original_neg != new_neg):
                #log(f"pH {ph:.2f} - {key}: (+{original_pos}, -{original_neg}) → (+{new_pos}, -{new_neg}) - pka: {pka_table[resname]}", silent)
                value[2] = new_pos
                value[3] = new_neg

    #log("\n", silent)
    return uncertainty_flags, local_contact_types


def resolve_uncertainty(name, uncertainty_flags, local_contact_types):
    flags = uncertainty_flags[name]
    new_atom_props = list(local_contact_types[name])    
    
    if flags.get('pos'):
        new_atom_props[2] = 1
    elif flags.get('neg'):
        new_atom_props[3] = 1

    return new_atom_props


def create_contacts(contacts_by_pair):
    count_types = {
        "hydrogen_bond":["HB",0],
        "hydrophobic":["HY",0],
        "attractive":["AT",0],
        "repulsive":["RE",0],
        "salt_bridge":["SB",0],
        "disulfide_bond":["DS",0],
        "stacking":["AS",0],
    }
    contacts = []

    for contact_list in contacts_by_pair.values():
        for entry in contact_list:
            
            contact_type = entry['type']
            is_uncertain = entry.get('is_uncertain', False)
            if is_uncertain:
                count_types[contact_type][1] += 1
                contact_type = f"uncertain_{contact_type}"
            elif contact_type.startswith('stacking'):
                count_types['stacking'][1] += 1
            else:
                count_types[contact_type][1] += 1
                
            contact = Contact(
                entry['protein_id'],
                entry['chain1'], entry['resnum1'], entry['resname1'], entry['atomname1'],
                entry['protein_id'],
                entry['chain2'], entry['resnum2'], entry['resname2'], entry['atomname2'],
                entry['distance'],
                contact_type,
                entry['atom1'],
                entry['atom2']
            )
            contacts.append(contact)

    return contacts, count_types
