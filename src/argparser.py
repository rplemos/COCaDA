"""
Author: Rafael Lemos - rafaellemos42@gmail.com
Date: 12/08/2024

License: MIT License
"""

from sys import exit
from argparse import ArgumentParser, ArgumentError, ArgumentTypeError
from multiprocessing import cpu_count
import re


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
        parser.add_argument('-d', '--distances', required=False, action='store_true', help='Processes custom contact distances based on the "contact_distances.txt" file.')

        args = parser.parse_args()

        files = args.files
        output = args.output
        distances = args.distances
                
        ncores = cpu_count()
        multi = args.multicore
        if multi is not None:
            if multi == 0:
                core = list(range(ncores))
            else:
                core = validate_core(multi, ncores)
        else:
            core = None
        
    except ArgumentError as e:
        print(f"Argument Error: {str(e)}")
        exit(1)

    except ValueError as e:
        print(f"Error: {str(e)}")
        exit(1)

    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}")
        exit(1)
    
    return files, core, output, distances
        
        
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
    # Check if it's a single core
    if value.isdigit():
        core = int(value)
        if core < 0 or core >= ncores:
            raise ArgumentTypeError(f"Core number {core} exceeds available cores (max: {ncores - 1})")
        return [core]
    
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