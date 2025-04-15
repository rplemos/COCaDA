"""
Author: Rafael Lemos - rafaellemos42@gmail.com
Date: 12/08/2024

License: MIT License
"""

from scripts.classes import Protein, Chain, Residue, Atom

import os
from numpy import mean, array
from numpy.linalg import svd


stacking = {
    'H':[10, 'CG','ND1','CE1','NE2','CD2'],
    'F':[11, 'CG','CD1','CE1','CZ','CE2','CD2'],
    'W':[14, 'CG','CD1','NE1','CE2','CZ2','CH2','CZ3','CE3','CD2'],
    'Y':[12, 'CG','CD1','CE1','CZ','CE2','CD2'],
}

residue_mapping = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}


def parse_pdb(pdb_file):
    """
    Parses a PDB file to create a Protein object.

    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        Protein: The Protein object populated with chain, residue, and atom objects.

    This function reads a PDB file line by line, extracts information about chains, residues,
    and atoms, and constructs a Protein object. It handles special cases such as alternative
    residue names, aromatic residues, and low-quality atoms. The function also assigns a unique
    identifier and title to the protein based on the header information.
    """
    
    current_protein = Protein()
    current_chain = None
    current_residue = None
    current_entity = None
    entity_chains = {}
    entity = None

    with open(pdb_file) as f:
        
        current_protein.id = os.path.basename(pdb_file).split(".")[0]
        
        for line in f:
            line = line.strip()
            
            if line == "ENDMDL":
                break
            
            # for interface checking
            elif line.startswith("COMPND"):
                if "MOL_ID" in line:
                    current_entity = line[-2]
                elif "CHAIN:" in line:
                    chains = line.split(":")[1].strip().replace(";","").replace(" ","")
                    entity_chains[current_entity] = chains.split(",")
        
            elif line.startswith("HEADER"):
                current_protein.id = line[62:]
                
            elif line.startswith("TITLE"):
                current_protein.set_title(line[10:])
                
            elif line.startswith("ATOM"):
                chain_id = line[21]
                
                if entity_chains:
                    for key, chains in entity_chains.items():
                        if chain_id in chains:
                            entity = key
                            break
                        
                resnum = int(line[22:26])
                # if resnum <= 0:
                #     continue
                resname = line[17:20]
                
                # alternative names for protonated histidines
                if resname in ["HID", "HIE", "HSP", "HSD", "HSE"]: 
                    resname = "HIS"
                
                if resname not in residue_mapping:
                    continue
                
                resname = residue_mapping.get(resname)                      

                if current_chain is None or current_chain.id != chain_id:  # new chain
                    residues = []
                    current_chain = Chain(chain_id, residues)
                    current_protein.chains.append(current_chain)
                    current_residue = None

                if current_residue is None:  # new residue
                    atoms = []
                    current_residue = Residue(resnum, resname, atoms, current_chain, False, None)
                    current_chain.residues.append(current_residue)
                
                if current_residue.resnum != resnum:
                    if len(current_residue.atoms) >= 1:
                        current_chain.residues.append(current_residue)
                    atoms = []
                    current_residue = Residue(resnum, resname, atoms, current_chain, False, None)
                                                                
                atomname = line[12:16].replace(" ", "")
                if atomname == "OXT" or atomname.startswith("H"): # OXT is the C-terminal Oxygen atom
                    continue
                
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                occupancy = float(line[55:60])
                
                if occupancy == 0 or occupancy >= 0.5: # ignores low quality atoms
                    if current_residue.atoms and current_residue.atoms[-1].atomname == atomname: # ignores the second one if both have occupancy == 0.5
                        continue
                    atom = Atom(atomname, x, y, z, occupancy, current_residue, entity) # creates atom
                    current_residue.atoms.append(atom)
                else:
                    continue

                # CHECKING FOR AROMATICS
                if current_residue.resname in stacking:
                    allowed = stacking[current_residue.resname][1:]
                    all_atoms_have_occupancy_one = all(atom.occupancy == 1 for atom in current_residue.atoms if atom.atomname in allowed)
                    
                    # if ring has only one conformation and the residue is complete (all atoms populated)
                    if all_atoms_have_occupancy_one and len(current_residue.atoms) == stacking[current_residue.resname][0]:
                        ring_atoms = array([[atom.x, atom.y, atom.z] for atom in current_residue.atoms if atom.atomname in stacking[current_residue.resname]])
                        
                        if ring_atoms.any():
                            centroid_atom = centroid(current_residue, ring_atoms, entity)
                            current_residue.atoms.append(centroid_atom)
                            current_residue.ring = True # flags the aromatic residue

                            normal_vector = calc_normal_vector(ring_atoms)
                            current_residue.normal_vector = normal_vector
                            
            elif line.startswith("END"):  
                # Handling cases where there is no ID
                if current_protein.id is None:
                    id = str(pdb_file).split("/")[-1]
                    id = id.split(".")[0]
                    current_protein.id = id  

    return current_protein


def parse_cif(cif_file):
    """
    Parses a .cif file to create a Protein object.

    Args:
        cif_file (str): Path to the cif file.

    Returns:
        Protein: The Protein object populated with chain, residue, and atom objects.

    This function reads a cif file, extracts information about chains, residues, and atoms,
    and constructs a Protein object. It handles specific fields for atom site data, processes
    only the first model for NMR files, and filters low-quality atoms. The function also assigns
    identifiers and titles based on the information from the file.
    """

    valid_atoms = ['N', 'C', 'O', 'S']
    
    current_protein = Protein()
    current_chain = None
    current_residue = None
    atomsite_block = False # _atom_site. lines
    atominfo_block = False # ATOM        lines
    atom_lines = []
    models = []
    title = None
    title_block = False

    with open(cif_file) as f:
        
        current_protein.id = os.path.basename(cif_file).split(".")[0]
        
        for line in f:
            line = line.strip()

            if line.startswith("_entry.id"):
                current_protein.id = line[-4:]

            # this title block can definitely be simplified, but there are a lot of edge cases to handle
            # current_protein.set_title() is set before returning the whole object
            if line.startswith("_struct.title"):
                title = line[len("_struct.title"):].strip()
                if title.startswith(";") and title.endswith(";"):
                    title = title[1:-1].strip()
                elif title.startswith("'") and title.endswith("'"):
                    title = title[1:-1].strip()
                elif title == "":
                    title_block = True
            elif title_block:
                if line.startswith(";"):
                    if title == "":
                        title = line[1:].strip()
                    else:
                        title_block = False
                else:
                    title += line.strip()
                    title_block = False

            if line.startswith("_atom_site.group_PDB"): # entering ATOM definition block
                atomsite_block = True
                line = line.split(".")[1]
                atom_lines.append(line)
                
            elif atomsite_block and line.startswith("_atom_site"):
                line = line.split(".")[1]
                atom_lines.append(line)
                
            elif atomsite_block and line.startswith("ATOM"): # maps the order of the columns                             
                atomname_index = atom_lines.index("label_atom_id")
                resname_index = atom_lines.index("label_comp_id")
                chain_index = atom_lines.index("label_asym_id")
                
                if "auth_seq_id" in atom_lines:
                    resnum_index = atom_lines.index("auth_seq_id")
                else:
                    resnum_index = atom_lines.index("label_seq_id")
                
                x_index = atom_lines.index("Cartn_x")
                y_index = atom_lines.index("Cartn_y")
                z_index = atom_lines.index("Cartn_z")
                occupancy_index = atom_lines.index("occupancy")
                model_index = atom_lines.index("pdbx_PDB_model_num")
                atom_element_index = atom_lines.index("type_symbol")
                entity_index = atom_lines.index("label_entity_id")
                                                               
                atomsite_block = False
                atominfo_block = True
                
            elif line.startswith("ATOM") and atominfo_block: # entering ATOM information block
                line = line.split()
                
                element = line[atom_element_index]
                if element not in valid_atoms:
                    continue
                
                models.append(int(line[model_index]))
                curr_model = int(line[model_index])
                if curr_model != models[0]: # parses only the first model (NMR files)
                    break
                    #return current_protein
                
                chain_id = line[chain_index]
                
                resnum = int(line[resnum_index])
                # if resnum <= 0:
                #     continue
                resname = line[resname_index]

                # alternative names for protonated histidines
                if resname in ["HID", "HIE", "HSP", "HSD", "HSE"]: 
                    resname = "HIS" 

                if resname not in residue_mapping:
                    continue
                
                resname = residue_mapping[resname]                            

                if current_chain is None or current_chain.id != chain_id:  # new chain
                    residues = []
                    current_chain = Chain(chain_id, residues)
                    current_protein.chains.append(current_chain)
                    current_residue = None

                if current_residue is None:  # first residue
                    atoms = []
                    current_residue = Residue(resnum, resname, atoms, current_chain, False, None)
                    current_chain.residues.append(current_residue)
                
                if current_residue.resnum != resnum: # new residue
                    if len(current_residue.atoms) >= 1:
                        current_chain.residues.append(current_residue) 
                    atoms = []
                    current_residue = Residue(resnum, resname, atoms, current_chain, False, None)
                                                                
                atomname = line[atomname_index]
                if atomname == "OXT"  or atomname.startswith("H"): # OXT is the C-terminal Oxygen atom
                    continue
                    
                x, y, z = float(line[x_index]), float(line[y_index]), float(line[z_index])
                occupancy = float(line[occupancy_index])
                
                entity = line[entity_index]
                    
                if (occupancy == 0 or occupancy >= 0.5): # ignores low quality atoms
                    if current_residue.atoms and current_residue.atoms[-1].atomname == atomname: # ignores the second one if both have occupancy == 0.5
                        continue
                    atom = Atom(atomname, x, y, z, occupancy, current_residue, entity) # creates atom
                    current_residue.atoms.append(atom)
                else:
                    continue
                                
                # CHECKING FOR AROMATICS
                if current_residue.resname in stacking:
                    allowed = stacking[current_residue.resname][1:]
                    all_atoms_have_occupancy_one = all(atom.occupancy == 1 for atom in current_residue.atoms if atom.atomname in allowed)
                    
                    # if ring has only one conformation and the residue is complete (all atoms populated)
                    if all_atoms_have_occupancy_one and len(current_residue.atoms) == stacking[current_residue.resname][0]:
                        ring_atoms = array([[atom.x, atom.y, atom.z] for atom in current_residue.atoms if atom.atomname in stacking[current_residue.resname]])
                        if ring_atoms.any():
                            centroid_atom = centroid(current_residue, ring_atoms, entity)
                            current_residue.atoms.append(centroid_atom)
                            current_residue.ring = True # flags the aromatic residue

                            normal_vector = calc_normal_vector(ring_atoms)
                            current_residue.normal_vector = normal_vector

            elif atominfo_block and line == "#":
                if resname in residue_mapping:
                    current_chain.residues.append(current_residue) # appends the last residue
                atominfo_block = False 
    
    if title is not None:
        current_protein.set_title(title.title().replace("'","").replace('"','').replace(",","."))
    else:
        current_protein.set_title(None)
    return current_protein


def centroid(residue, ring_atoms, entity):
    """
    Calculates the centroid of a set of ring atoms and creates a centroid Atom.

    Args:
        residue (Residue): The residue containing the ring atoms.
        ring_atoms (array): Array of ring atom coordinates.

    Returns:
        Atom: A new Atom object representing the centroid of the ring.

    The function computes the centroid of the given ring atoms and creates a new Atom object
    with this centroid position. The atom is labeled as "RNG" for ring centroid.
    """
    
    centroid = mean(ring_atoms, axis = 0)
    centroid_atom = Atom("RNG", centroid[0], centroid[1], centroid[2], 1, residue, entity)
    
    return centroid_atom


def calc_normal_vector(ring_atoms):
    """
    Computes the normal vector to the plane of a set of ring atoms.

    Args:
        ring_atoms (array): Array of ring atom coordinates.

    Returns:
        array: The normal vector to the plane of the ring atoms.

    The function calculates the normal vector to the plane defined by the given ring atoms
    using Singular Value Decomposition (SVD). The normal vector is extracted from the last
    row of the V^T matrix obtained from SVD.
    """
    
    centroid = mean(ring_atoms, axis = 0) # axis=0 -> mean through the columns
    centered_ring_atoms = ring_atoms - centroid # normalizes to origin
    
    # Use SVD to calculate the plane
    _, _, vh = svd(centered_ring_atoms) # vh = V^T
    normal_vector = vh[2]  # normal vector is the last row of the V^T matrix
    
    return normal_vector
