"""
Author: Rafael Lemos - rafaellemos42@gmail.com
Date: 12/08/2024

License: MIT License
"""

from sys import exit
from argparse import ArgumentParser, ArgumentError, ArgumentTypeError
from multiprocessing import cpu_count
import re
import json


def cl_parse():
    """
    Parses command-line arguments for a PDB/mmCIF parser and contact detection tool.

    Returns:
        tuple: A tuple containing the parsed values:
            - files (list): List of input files.
            - multicore (bool): Select MultiCore mode.
            - core (int): Select cores to use.
            - output (bool): Whether to output results to files.

    Raises:
        ArgumentError: If there's an issue with the command-line arguments.
        ValueError: If an invalid processing mode is specified.
        Exception: For any other unexpected errors during argument parsing.
    """
    
    try:
        parser = ArgumentParser(description='COCαDA - Large-Scale Protein Interatomic Contact Cutoff Optimization by Cα Distance Matrices.')
        parser.add_argument('-f', '--files', nargs='+', required=True, type=validate_file, help='List of files in pdb/cif format (at least one required). Wildcards are accepted (ex. -f *.cif).')
        parser.add_argument('-m', '--multicore', required=False, nargs='?', const=0, help='Use MultiCore mode. Default uses all available cores, and selections can be defined based on the following: -m X = specific single core. -m X-Y = range of cores from X to Y. -m X,Y,Z... = specific multiple cores.')
        parser.add_argument('-o', '--output', required=False, nargs='?', const='./outputs', help='Outputs the results to files in the given folder. Default is ./outputs.')
        parser.add_argument('-r', '--region', required=False, nargs='?', help='Define only a region of residues to be analyzed. Selections can be defined based on the following: -r X-Y = range of residues from X to Y. -r X,Y,Z... = specific multiple residues.')
        parser.add_argument('-i', '--interface', required=False, nargs='?', const='interface.csv', help='Calculate only interface contacts.')        
        parser.add_argument('-d', '--distances', nargs='?', type=validate_distances, default=False, help='Processes custom contact distances based on the "contact_distances.txt" file.')
        parser.add_argument('-ph', '--ph', type=validate_ph, default=None, help='pH value (0-14)')
        parser.add_argument('-s', '--silent', required=False, action='store_true', help='Suppresses non-essential console output.')
        parser.add_argument('-c', '--chains', required=False, nargs='?', help='Define only specific chains to be analyzed. -c A = only one chain. -c A,B,C... = specific multiple chains.')

        args = parser.parse_args()

        files = args.files
        output = args.output
        interface = args.interface
        distances = args.distances
        ph = args.ph
        silent = args.silent
                
        ncores = cpu_count()
        multi = args.multicore
        if multi is not None:
            if multi == 0:
                core = list(range(ncores))
            else:
                core = validate_core(multi, ncores)
        else:
            core = None
            
        region_values = args.region
        if region_values is not None:
            region = validate_region(region_values)
        else:
            region = None
            
        chain_values = args.chains
        if chain_values is not None:
            chains = validate_chains(chain_values)
        else:
            chains = None

    except ArgumentError as e:
        print(f"Argument Error: {str(e)}")
        exit(1)

    except ValueError as e:
        print(f"Error: {str(e)}")
        exit(1)

    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}")
        exit(1)
    
    return files, core, output, region, chains, interface, distances, ph, silent
        
        
def validate_file(value):
    """
    Validates a file path to ensure it has a proper extension for PDB or mmCIF files.

    If the file has a valid extension, the function returns the file path. Otherwise, it raises an `ArgumentTypeError`.

    Args:
        value (str): The file path to validate.

    Returns:
        str: The validated file path.

    Raises:
        ArgumentTypeError: If the file does not have a valid extension.
    """
    
    if value.endswith('.pdb') or value.endswith('.cif'):
        return value
    else:
        raise ArgumentTypeError(f"{value} is not a valid file. File must end with '.pdb' or '.cif'")



def validate_core(value, ncores):
    """
    Validates the --core argument to ensure it follows the correct format.
    Supports single core, range of cores, and list of cores.

    Args:
        value (str): The value input by the user for the --core argument.
        ncores (int): The maximum number of cores on the system.

    Returns:
        list: A list of valid cores to use.

    Raises:
        ArgumentTypeError: If the input is not valid or exceeds available cores.
    """
    # Check if it's a single number representing the number of cores to use
    if value.isdigit():
        core_count = int(value)
        if core_count <= 0 or core_count > ncores:
            raise ArgumentTypeError(f"Requested number of cores {core_count} exceeds available cores (max: {ncores})")
        return core_count
    
    # Check if it's a range (e.g. 10-19)
    range_match = re.match(r'^(\d+)-(\d+)$', value)
    if range_match:
        start_core, end_core = map(int, range_match.groups())
        if start_core < 0 or end_core >= ncores or start_core > end_core:
            raise ArgumentTypeError(f"Invalid range {start_core}-{end_core}, ensure it's within [0-{ncores - 1}]")
        return list(range(start_core, end_core + 1))

    # Check if it's a list of cores (e.g. 10,32,65)
    list_match = re.match(r'^(\d+(,\d+)+)$', value)
    if list_match:
        core_list = list(map(int, value.split(',')))
        if any(core < 0 or core >= ncores for core in core_list):
            raise ArgumentTypeError(f"One or more cores exceed available cores (max: {ncores - 1})")
        return core_list
    
    raise ArgumentTypeError(f"Invalid core format: {value}. Use a single core, a range (x-y), or a list (x,y,z).")


def validate_region(region):
    """
    Validates and parses a region input, which can be either a range (e.g., "10-19")
    or a comma-separated list of values (e.g., "10,32,65").

    Args:
        region (str): The input string representing the region selection.

    Returns:
        list: A list of integers representing the parsed region values.

    Raises:
        ArgumentTypeError: If the input format is invalid or contains negative values.
    """
    
    # Check if it's a range (e.g. 10-19)
    range_match = re.match(r'^(\d+)-(\d+)$', region)
    if range_match:
        start_res, end_res = map(int, range_match.groups())
        if start_res < 0 or start_res > end_res:
            raise ArgumentTypeError(f"Invalid range {start_res}-{end_res}]")
        return list(range(start_res, end_res + 1))

    # Check if it's a list (e.g. 10,32,65)
    list_match = re.match(r'^(\d+(,\d+)+)$', region)
    if list_match:
        res_list = list(map(int, region.split(',')))
        if any(core < 0 for core in res_list):
            raise ArgumentTypeError("One or more value is not valid.")
        return res_list
    
    # Check if it's a list of single uppercase letters (e.g. A,B,C)
    letter_list_match = re.match(r'^([a-zA-Z](,[a-zA-Z])*)$', region)
    if letter_list_match:
        return region.split(',')
    
    raise ArgumentTypeError(f"Invalid region format: {region}. Use a range (x-y) or a list (x,y,z).")


def validate_chains(chains):
    # Check if it's a list of single uppercase letters (e.g. A,B,C)
    letter_list_match = re.match(r'^([a-zA-Z](,[a-zA-Z])*)$', chains)
    if letter_list_match:
        return chains.split(',')
    
    raise ArgumentTypeError(f"Invalid region format: {chains}. Use one value (A) or a list (A,B,C...).")


def validate_ph(value):
    """
    Validates the pH argument from the command line, ensuring that the provided value 
    is a valid floating-point number within the acceptable range (0.0 to 14.0).

    Parameters:
        value (str): The pH value as a string (from command-line input).

    Returns:
        float: The validated pH value as a float.

    Raises:
        ArgumentTypeError: If the input is not a number or is outside the 0–14 range.
    """
    
    try:
        ph = float(value)
    except ValueError:
        raise ArgumentTypeError(f"Invalid pH value: '{value}' is not a number.")
    if not (0.0 <= ph <= 14.0):
        raise ArgumentTypeError(f"Invalid pH: Must be between 0 and 14 (got {ph}).")
    return ph


def validate_distances(value):
    distance_keys = [
        'salt_bridge',
        'hydrophobic',
        'hydrogen_bond',
        'repulsive',
        'attractive',
        'disulfide_bond',
        'aromatic'
    ]

    print(value)
    print(type(value))
    loaded_distances = {}
    
    def validate_categories(categories): 
        for key, (min_val, max_val) in categories.items():
            if min_val < 0 or max_val < 0:
                raise ValueError(f"Invalid values for '{key}': values must be positive.")
            if min_val >= max_val:
                raise ValueError(f"Invalid range for '{key}': min ({min_val}) must be less than max ({max_val}).")
        return categories

    if value is True:  # User passed just '-d' with no value → JSON mode
        try:
            with open("./contact_distances.json", "r") as f:
                loaded_distances = json.load(f)
        except Exception as e:
            raise ArgumentTypeError(f"Error loading JSON file: {e}")
    else:  # User passed 14 comma-separated values
        parts = value.split(',')
        if len(parts) != 14:
            raise ArgumentTypeError(f"Expected 14 comma-separated distance values, got {len(parts)}.")
        try:
            floats = [float(x) for x in parts]
        except ValueError:
            raise ArgumentTypeError(f"All distance values must be numeric floats. Got: {parts}")
        
        for i, key in enumerate(distance_keys):
            start = i*2
            loaded_distances[key] = [floats[start], floats[start + 1]]
     
    try:
        validated_distances = validate_categories({key: tuple(value) for key, value in loaded_distances.items()})  
    except ValueError:
        raise ArgumentTypeError("Invalid!")

    return validated_distances