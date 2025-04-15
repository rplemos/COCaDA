"""
Author: Rafael Lemos - rafaellemos42@gmail.com
Date: 12/08/2024

License: MIT License
"""

from math import dist
from numpy import dot, arccos, degrees
from numpy.linalg import norm

from scripts.classes import Contact
from scripts.distances import distances
import scripts.conditions as conditions


def contact_detection(protein, region, interface, custom_distances, epsilon):
    """
    Detects contacts between atoms in a given protein.

    Args:
        protein (Protein): The protein object containing chain, residue and atom objects.

    Returns:
        list: A list of Contact objects representing the detected contacts.
    """

    residues = list(protein.get_residues())
    contacts = []
    interface_res = set()
    max_ca_distance = 20.47 # 0.01 higher than the Arg-Arg pair
    
    categories = custom_distances if custom_distances else conditions.categories
    if epsilon > 0:
        max_ca_distance += epsilon
        updated_distances = {key: value + epsilon for key, value in distances.items()}
    else:
        updated_distances = distances
        
    for i, residue1 in enumerate(residues[1:]):
        for _, residue2 in enumerate(residues[i+1:], start=i+1):
            
            if residue1.resnum == residue2.resnum and residue1.chain.id == residue2.chain.id: # ignores same residue
                continue
            
            if region and (residue1.resnum not in region or residue2.resnum not in region):
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

                    contact = Contact(protein.id, residue1.chain.id, residue1.resnum, residue1.resname, ring1.atomname, 
                                    protein.id, residue2.chain.id, residue2.resnum, residue2.resname, ring2.atomname, 
                                    float(f"{distance:.2f}"), "stacking"+stack_type, ring1, ring2)
                    
                    contacts.append(contact)
                    
            for atom1 in residue1.atoms:
                for atom2 in residue2.atoms:
                    
                    if interface and atom1.entity == atom2.entity:
                        continue
                    
                    name1 = f"{atom1.residue.resname}:{atom1.atomname}" # matches the pattern from conditions dictionary
                    name2 = f"{atom2.residue.resname}:{atom2.atomname}"

                    if name1 in conditions.contact_types and name2 in conditions.contact_types: # excludes the RNG atom and any different other
                        
                        distance = dist((atom1.x, atom1.y, atom1.z), (atom2.x, atom2.y, atom2.z))
                        
                        if distance <= 6: # max distance for contacts

                            for contact_type, distance_range in categories.items():

                                if contact_type == 'hydrogen_bond' and (abs(residue2.resnum - residue1.resnum) <= 3): # skips alpha-helix for h-bonds
                                    continue
                                
                                if distance_range[0] <= distance <= distance_range[1]: # fits the range
                                    if conditions.contact_conditions[contact_type](name1, name2): # fits the type of contact
                                                                                                
                                        contact = Contact(protein.id, residue1.chain.id, residue1.resnum, residue1.resname, atom1.atomname, 
                                                        protein.id, residue2.chain.id, residue2.resnum, residue2.resname, atom2.atomname, 
                                                        float(f"{distance:.2f}"), contact_type, atom1, atom2)

                                        contacts.append(contact)
                                        
                                        interface_res.add(f"{residue1.chain.id},{residue1.resnum},{residue1.resname}")
    return contacts, interface_res


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
