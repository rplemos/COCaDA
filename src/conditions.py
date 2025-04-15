"""
Author: Rafael Lemos - rafaellemos42@gmail.com
Date: 12/08/2024

License: MIT License
"""

# RULES
# 1 - must be made by different residue atoms
# 2 - aromatic = aromatic + aromatic
# 3 - hydrogenb => aceptor + donor
# 4 - hydrophobic: hydrofobic + hydrofobic
# 5 - Repulsive: positive=>positive or negative=>negative
# 6 - Atractive: positive=>negative or negative=>positive
# 7 - salt_bridge: positive=>negative or negative=>positive

categories = {
    'salt_bridge': (0, 3.9), # minimum and maximum distances, in Angstroms
    'hydrophobic': (2, 4.5),
    'hydrogen_bond': (0, 3.9),
    'repulsive': (2, 6),
    'attractive': (3.9, 6),
    'disulfide_bond': (0, 2.8),
    'aromatic': (2, 5)
}


contact_conditions = {
    'salt_bridge': lambda name1, name2: (contact_types[name1][2] == 1 and contact_types[name2][3] == 1) or (contact_types[name1][3] == 1 and contact_types[name2][2] == 1),
    'disulfide_bond': lambda name1, name2: name1 == "C:SG" and name2 == "C:SG",
    'hydrogen_bond': lambda name1, name2: ((contact_types[name1][4] == 1 and contact_types[name2][5] == 1) or (contact_types[name1][5] == 1 and contact_types[name2][4] == 1)),  
    'hydrophobic': lambda name1, name2: contact_types[name1][0] == 1 and contact_types[name2][0] == 1,
    'repulsive': lambda name1, name2: (contact_types[name1][2] == 1 and contact_types[name2][2] == 1) or (contact_types[name1][3] == 1 and contact_types[name2][3] == 1),
    'attractive': lambda name1, name2: (contact_types[name1][2] == 1 and contact_types[name2][3] == 1) or (contact_types[name1][3] == 1 and contact_types[name2][2] == 1),
    'aromatic': lambda name1, name2: (contact_types[name1][1] == 2 and contact_types[name2][1] == 2) # aromatics are checked elsewhere
}


# 'RES:ATOM':	[	Hydrophobic,	Aromatic,	Positive,	Negative,	Donor,	Acceptor	]
# 'A:CA':		[	0|1,			0|1,		0|1,		0|1,		0|1,	0|1			]
contact_types = { 
    'A:N':[0, 0, 0, 0, 1, 0],
    'A:CA':[0, 0, 0, 0, 0, 0],
    'A:C':[0, 0, 0, 0, 0, 0],
    'A:O':[0, 0, 0, 0, 0, 1],
    'A:CB':[1, 0, 0, 0, 0, 0],
    
    'R:N':[0, 0, 0, 0, 1, 0],
    'R:CA':[0, 0, 0, 0, 0, 0],
    'R:C':[0, 0, 0, 0, 0, 0],
    'R:O':[0, 0, 0, 0, 0, 1],
    'R:CB':[1, 0, 0, 0, 0, 0],
    'R:CG':[1, 0, 0, 0, 0, 0],
    'R:CD':[0, 0, 0, 0, 0, 0],
    'R:NE':[0, 0, 1, 0, 1, 0],
    'R:CZ':[0, 0, 1, 0, 0, 0],
    'R:NH1':[0, 0, 1, 0, 1, 0],
    'R:NH2':[0, 0, 1, 0, 1, 0],
    
    'N:N':[0, 0, 0, 0, 1, 0],
    'N:CA':[0, 0, 0, 0, 0, 0],
    'N:C':[0, 0, 0, 0, 0, 0],
    'N:O':[0, 0, 0, 0, 0, 1],
    'N:CB':[1, 0, 0, 0, 0, 0],
    'N:CG':[0, 0, 0, 0, 0, 0],
    'N:OD1':[0, 0, 0, 0, 0, 1],
    'N:ND2':[0, 0, 0, 0, 1, 0],
    
    'D:N':[0, 0, 0, 0, 1, 0],
    'D:CA':[0, 0, 0, 0, 0, 0],
    'D:C':[0, 0, 0, 0, 0, 0],
    'D:O':[0, 0, 0, 0, 0, 1],
    'D:CB':[1, 0, 0, 0, 0, 0],
    'D:CG':[0, 0, 0, 0, 0, 0],
    'D:OD1':[0, 0, 0, 1, 0, 1],
    'D:OD2':[0, 0, 0, 1, 0, 1],
    
    'C:N':[0, 0, 0, 0, 1, 0],
    'C:CA':[0, 0, 0, 0, 0, 0],
    'C:C':[0, 0, 0, 0, 0, 0],
    'C:O':[0, 0, 0, 0, 0, 1],
    'C:CB':[1, 0, 0, 0, 0, 0],
    'C:SG':[0, 0, 0, 0, 1, 1],
    
    'Q:N':[0, 0, 0, 0, 1, 0],
    'Q:CA':[0, 0, 0, 0, 0, 0],
    'Q:C':[0, 0, 0, 0, 0, 0],
    'Q:O':[0, 0, 0, 0, 0, 1],
    'Q:CB':[1, 0, 0, 0, 0, 0],
    'Q:CG':[1, 0, 0, 0, 0, 0],
    'Q:CD':[0, 0, 0, 0, 0, 0],
    'Q:OE1':[0, 0, 0, 0, 0, 1],
    'Q:NE2':[0, 0, 0, 0, 1, 0],
    
    'E:N':[0, 0, 0, 0, 1, 0],
    'E:CA':[0, 0, 0, 0, 0, 0],
    'E:C':[0, 0, 0, 0, 0, 0],
    'E:O':[0, 0, 0, 0, 0, 1],
    'E:CB':[1, 0, 0, 0, 0, 0],
    'E:CG':[1, 0, 0, 0, 0, 0],
    'E:CD':[0, 0, 0, 0, 0, 0],
    'E:OE1':[0, 0, 0, 1, 0, 1],
    'E:OE2':[0, 0, 0, 1, 0, 1],
    
    'G:N':[0, 0, 0, 0, 1, 0],
    'G:CA':[0, 0, 0, 0, 0, 0],
    'G:C':[0, 0, 0, 0, 0, 0],
    'G:O':[0, 0, 0, 0, 0, 1],
    
    'H:N':[0, 0, 0, 0, 1, 0],
    'H:CA':[0, 0, 0, 0, 0, 0],
    'H:C':[0, 0, 0, 0, 0, 0],
    'H:O':[0, 0, 0, 0, 0, 1],
    'H:CB':[1, 0, 0, 0, 0, 0],
    'H:CG':[0, 1, 0, 0, 0, 0],
    'H:ND1':[0, 1, 1, 0, 1, 1],
    'H:CD2':[0, 1, 0, 0, 0, 0],
    'H:CE1':[0, 1, 0, 0, 0, 0],
    'H:NE2':[0, 1, 1, 0, 1, 1],
    
    'I:N':[0, 0, 0, 0, 1, 0],
    'I:CA':[0, 0, 0, 0, 0, 0],
    'I:C':[0, 0, 0, 0, 0, 0],
    'I:O':[0, 0, 0, 0, 0, 1],
    'I:CB':[1, 0, 0, 0, 0, 0],
    'I:CG1':[1, 0, 0, 0, 0, 0],
    'I:CG2':[1, 0, 0, 0, 0, 0],
    'I:CD1':[1, 0, 0, 0, 0, 0],
    
    'L:N':[0, 0, 0, 0, 1, 0],
    'L:CA':[0, 0, 0, 0, 0, 0],
    'L:C':[0, 0, 0, 0, 0, 0],
    'L:O':[0, 0, 0, 0, 0, 1],
    'L:CB':[1, 0, 0, 0, 0, 0],
    'L:CG':[1, 0, 0, 0, 0, 0],
    'L:CD1':[1, 0, 0, 0, 0, 0],
    'L:CD2':[1, 0, 0, 0, 0, 0],
    
    'K:N':[0, 0, 0, 0, 1, 0],
    'K:CA':[0, 0, 0, 0, 0, 0],
    'K:C':[0, 0, 0, 0, 0, 0],
    'K:O':[0, 0, 0, 0, 0, 1],
    'K:CB':[1, 0, 0, 0, 0, 0],
    'K:CG':[1, 0, 0, 0, 0, 0],
    'K:CD':[1, 0, 0, 0, 0, 0],
    'K:CE':[0, 0, 0, 0, 0, 0],
    'K:NZ':[0, 0, 1, 0, 1, 0],
    
    'M:N':[0, 0, 0, 0, 1, 0],
    'M:CA':[0, 0, 0, 0, 0, 0],
    'M:C':[0, 0, 0, 0, 0, 0],
    'M:O':[0, 0, 0, 0, 0, 1],
    'M:CB':[1, 0, 0, 0, 0, 0],
    'M:CG':[1, 0, 0, 0, 0, 0],
    'M:SD':[0, 0, 0, 0, 0, 1],
    'M:CE':[1, 0, 0, 0, 0, 0],
    
    'F:N':[0, 0, 0, 0, 1, 0],
    'F:CA':[0, 0, 0, 0, 0, 0],
    'F:C':[0, 0, 0, 0, 0, 0],
    'F:O':[0, 0, 0, 0, 0, 1],
    'F:CB':[1, 0, 0, 0, 0, 0],
    'F:CG':[1, 1, 0, 0, 0, 0],
    'F:CD1':[1, 1, 0, 0, 0, 0],
    'F:CD2':[1, 1, 0, 0, 0, 0],
    'F:CE1':[1, 1, 0, 0, 0, 0],
    'F:CE2':[1, 1, 0, 0, 0, 0],
    'F:CZ':[1, 1, 0, 0, 0, 0],
    
    'P:N':[0, 0, 0, 0, 0, 0],
    'P:CA':[0, 0, 0, 0, 0, 0],
    'P:C':[0, 0, 0, 0, 0, 0],
    'P:O':[0, 0, 0, 0, 0, 1],
    'P:CB':[1, 0, 0, 0, 0, 0],
    'P:CG':[1, 0, 0, 0, 0, 0],
    'P:CD':[0, 0, 0, 0, 0, 0],
    
    'S:N':[0, 0, 0, 0, 1, 0],
    'S:CA':[0, 0, 0, 0, 0, 0],
    'S:C':[0, 0, 0, 0, 0, 0],
    'S:O':[0, 0, 0, 0, 0, 1],
    'S:CB':[0, 0, 0, 0, 0, 0],
    'S:OG':[0, 0, 0, 0, 1, 1],
    
    'T:N':[0, 0, 0, 0, 1, 0],
    'T:CA':[0, 0, 0, 0, 0, 0],
    'T:C':[0, 0, 0, 0, 0, 0],
    'T:O':[0, 0, 0, 0, 0, 1],
    'T:CB':[0, 0, 0, 0, 0, 0],
    'T:OG1':[0, 0, 0, 0, 1, 1],
    'T:CG2':[1, 0, 0, 0, 0, 0],
    
    'W:N':[0, 0, 0, 0, 1, 0],
    'W:CA':[0, 0, 0, 0, 0, 0],
    'W:C':[0, 0, 0, 0, 0, 0],
    'W:O':[0, 0, 0, 0, 0, 1],
    'W:CB':[1, 0, 0, 0, 0, 0],
    'W:CG':[1, 1, 0, 0, 0, 0],
    'W:CD1':[0, 1, 0, 0, 0, 0],
    'W:CD2':[1, 1, 0, 0, 0, 0],
    'W:NE1':[0, 1, 0, 0, 1, 0],
    'W:CE2':[0, 1, 0, 0, 0, 0],
    'W:CE3':[1, 1, 0, 0, 0, 0],
    'W:CZ2':[1, 1, 0, 0, 0, 0],
    'W:CZ3':[1, 1, 0, 0, 0, 0],
    'W:CH2':[1, 1, 0, 0, 0, 0],
    
    'Y:N':[0, 0, 0, 0, 1, 0],
    'Y:CA':[0, 0, 0, 0, 0, 0],
    'Y:C':[0, 0, 0, 0, 0, 0],
    'Y:O':[0, 0, 0, 0, 0, 1],
    'Y:CB':[1, 0, 0, 0, 0, 0],
    'Y:CG':[1, 1, 0, 0, 0, 0],
    'Y:CD1':[1, 1, 0, 0, 0, 0],
    'Y:CD2':[1, 1, 0, 0, 0, 0],
    'Y:CE1':[1, 1, 0, 0, 0, 0],
    'Y:CE2':[1, 1, 0, 0, 0, 0],
    'Y:CZ':[0, 1, 0, 0, 0, 0],
    'Y:OH':[0, 0, 0, 0, 1, 1],
    
    'V:N':[0, 0, 0, 0, 1, 0],
    'V:CA':[0, 0, 0, 0, 0, 0],
    'V:C':[0, 0, 0, 0, 0, 0],
    'V:O':[0, 0, 0, 0, 0, 1],
    'V:CB':[1, 0, 0, 0, 0, 0],
    'V:CG1':[1, 0, 0, 0, 0, 0],
    'V:CG2':[1, 0, 0, 0, 0, 0]
}
