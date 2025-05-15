"""
Author: Rafael Lemos - rafaellemos42@gmail.com
Date: 12/08/2024

License: MIT License
"""

import os
import json
from timeit import default_timer as timer

import src.argparser as argparser
import src.classes as classes
import src.process as process


def main():
    """
    Main function for the script.

    This function parses command-line arguments, sets up the environment based on the specified mode,
    and runs the appropriate processing function (single-core or multi-core) on the input files.
    It also manages core affinity, output folder creation, and timing of the entire process.
    """
    global_time_start = timer()
    
    file_list, core, output, custom_distances = argparser.cl_parse()
    
    print("--------------COCaDA----------------\n")
    
    # context object for shared parameters
    context = classes.ProcessingContext(core=core, output=output, custom_distances=custom_distances) 
    
    if core is not None:  # Set specific core affinity
        print("Multicore mode selected.")
    else:
        print("Running on single mode with no specific core.")
               
    if output:
        print(f"Generating outputs in '{output}' folder.")
        if not os.path.exists(output):
            os.makedirs(output)
    else:
        output = None
        
    if custom_distances:
        print("Using custom distances provided by the user.")
        with open("./contact_distances.json","r") as f:
            loaded_distances = json.load(f)
        try:
            validated_distances = process.validate_categories({key: tuple(value) for key, value in loaded_distances.items()})
            max_value = max(y for x in validated_distances.values() for y in x)
            if max_value > 6:
                context.epsilon = max_value - 6
        except ValueError as e:
            print(e)  
            exit(1)
            
        context.custom_distances = validated_distances

    print()
    process_func = process.single if core is None else process.multi_batch
    process_func(file_list, context)
    
    print("\n------------------------------------\n")
    print(f"Total time elapsed: {(timer() - global_time_start):.3f}s\n")


if __name__ == "__main__":
    main()
