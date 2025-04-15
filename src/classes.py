"""
Author: Rafael Lemos - rafaellemos42@gmail.com
Date: 12/08/2024

License: MIT License
"""

class ProcessingContext:
    """
    Stores processing parameters passed by command-line.

    Attributes:
        core (str or None): The core selection for processing. Defaults to None.
        output (str or None): The output file path or identifier. Defaults to None.
    """

    def __init__(self, core=None, output=None, custom_distances=False, epsilon=0):
        self.core = core
        self.output = output
        self.custom_distances = custom_distances
        self.epsilon = epsilon

class Protein:
    """
    Represents a protein structure, including its title, ID, and chains.

    Attributes:
        title (str or None): The title of the protein.
        id (str or None): The ID of the protein.
        chains (list of Chain): A list of Chain objects representing the chains in the protein.
    """
    
    
    def __init__(self):
        """
        Initializes a new Protein instance with default values.
        """
        
        self.title = None
        self.id = None
        self.chains = []

    def set_title(self, title):
        """
        Sets the title of the protein. If a title is already set, appends the new title.

        Args:
            title (str): The title to set or append to the protein's title.
        """
        
        if self.title is None:
            self.title = title
        else:
            self.title += " " + title.strip()

    def get_chains(self):
        """
        Generator that yields each chain in the protein.

        Yields:
            Chain: A chain object in the protein.
        """
        
        for chain in self.chains:
            yield chain

    def get_residues(self):
        """
        Generator that yields each residue in all chains of the protein.

        Yields:
            Residue: A residue object in the protein's chains.
        """
        
        for chain in self.get_chains():
            for residue in chain.residues:
                yield residue
                
    def true_count(self):
        """
        Counts the total number of residues in the protein.

        Returns:
            int: The total number of residues in the protein.
        """
        
        return sum(1 for _ in self.get_residues())
            

class Chain:
    """
    Represents a chain in a protein, containing residues.

    Attributes:
        id (str): The ID of the chain.
        residues (list of Residue): A list of Residue objects representing the residues in the chain.
    """

    def __init__(self, id, residues):
        """
        Initializes a new Chain instance.
        """
        
        self.id = id
        self.residues = residues


class Residue:
    """
    Represents a residue in a protein chain.

    Attributes:
        resnum (int): The residue number.
        resname (str): The residue name (e.g., 'ALA' for alanine).
        atoms (list of Atom): A list of Atom objects in the residue.
        chain (Chain): The chain to which the residue belongs.
        ring (bool): Indicates whether the residue has a ring structure.
        normal_vector (tuple): The normal vector associated with the residue.
    """
    
    def __init__(self, resnum, resname, atoms, chain, ring, normal_vector):
        """
        Initializes a new Residue instance.
        """
        
        self.resnum = resnum
        self.resname = resname
        self.atoms = atoms
        self.chain = chain
        self.ring = ring
        self.normal_vector = normal_vector
  
      
class Atom:
    """
    Represents an atom in a residue.

    Attributes:
        atomname (str): The name of the atom (e.g., 'CA' for alpha carbon).
        x (float): The x-coordinate of the atom.
        y (float): The y-coordinate of the atom.
        z (float): The z-coordinate of the atom.
        occupancy (float): The occupancy value of the atom.
        residue (Residue): The residue to which the atom belongs.
        entity (int): The entity in which the atom is located.
    """
    
    def __init__(self, atomname, x, y, z, occupancy, residue, entity):
        """
        Initializes a new Atom instance.
        """
        
        self.atomname = atomname
        self.x = x
        self.y = y
        self.z = z
        self.occupancy = occupancy
        self.residue = residue
        self.entity = entity


class Contact: 
    """
    Represents a contact between two atoms in different residues.

    Attributes:
        id1 (str): The ID of the first residue.
        chain1 (str): The chain of the first residue.
        residue_num1 (int): The residue number of the first residue.
        residue_name1 (str): The residue name of the first residue.
        atom1 (str): The atom name of the first atom.
        id2 (str): The ID of the second residue.
        chain2 (str): The chain of the second residue.
        residue_num2 (int): The residue number of the second residue.
        residue_name2 (str): The residue name of the second residue.
        atom2 (str): The atom name of the second atom.
        distance (float): The distance between the two atoms.
        type (str): The type of contact (e.g., hydrogen bond, hydrophobic).
        atom_object1 (Atom): The Atom object representing the first atom.
        atom_object2 (Atom): The Atom object representing the second atom.
    """
    
    def __init__(self, id1, chain1, residue_num1, residue_name1, atom1, 
                 id2, chain2, residue_num2, residue_name2, atom2, 
                 distance, type, atom_object1, atom_object2):
        """
        Initializes a new Contact instance.
        """
        
        self.id1 = id1
        self.chain1 = chain1
        self.residue_num1 = residue_num1
        self.residue_name1 = residue_name1
        self.atom1 = atom1
        self.id2 = id2
        self.chain2 = chain2
        self.residue_num2 = residue_num2
        self.residue_name2 = residue_name2
        self.atom2 = atom2
        self.distance = distance
        self.type = type
        self.atom_object1 = atom_object1
        self.atom_object2 = atom_object2
    
    def print_text(self):
        """
        Generates a formatted string describing the contact between two atoms.

        Returns:
            str: A string describing the contact, including chain, residue, atom information, and distance.
        """
        
        map_type = {
            "hydrogen_bond":"HB",
            "hydrophobic":"HY",
            "attractive":"AT",
            "repulsive":"RE",
            "salt_bridge":"SB",
            "disulfide_bond":"DS",
            "stacking-other":"AS",
            "stacking-parallel":"AS", # on v.1 all aromatic stackings will be considered the same
            "stacking-perpendicular":"AS" # need to reimplement later
        }
        
        all_values = list(self.__dict__.values())
        all_values[11] = map_type[all_values[11]]
        return f"{all_values[1]},{all_values[2]},{all_values[3]},{all_values[4]},{all_values[6]},{all_values[7]},{all_values[8]},{all_values[9]},{all_values[10]},{all_values[11]}"
    
        # verbose
        # return f"{all_values[1]}-{all_values[2]}{all_values[3]}:{all_values[4]} and {all_values[6]}-{all_values[7]}{all_values[8]}:{all_values[9]}: {all_values[10]} A. {all_values[11].capitalize()}"

