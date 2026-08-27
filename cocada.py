"""
COCaDA v1.6
Author: Rafael Lemos - rafaellemos@ufmg.br
Date: 26/03/2026

License: MIT License
"""

import os
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
    
    file_list, core, output, region, chains, interface, distances, modified_distances, ph, silent, interchain, web, propedia = argparser.cl_parse()
    
    process.log("\n--------------COCaDA----------------\n", silent)
    
    # context object for shared parameters
    context = classes.ProcessingContext(core=core, output=output, region=region, chains=chains, interface=interface, 
                                        distances=distances, modified_distances=modified_distances, ph=ph, silent=silent, interchain=interchain, web=int(web), propedia=propedia)

    if core is not None:  # Set specific core affinity
        process.log("Multicore mode selected.", silent)
    else:
        process.log("Running on single mode with no specific core.", silent)

    if interface:
        process.log("Calculating only interface contacts.", silent)
        
    if region:
        process.log(f"Calculating contacts in the region: {region}", silent)
    
    if chains:
        process.log(f"Calculating contacts in the chains: {chains}", silent)
        
    if interchain:
        process.log(f"Calculating only interchain contacts.", silent)
               
    if output:
        process.log(f"Generating outputs in '{output}' folder.", silent)
        if not os.path.exists(output):
            os.makedirs(output)
    else:
        output = None
        
    if ph and ph != -1:
        process.log(f"Changing protonation states of pH-sensitive atoms using pH value of {ph}.\n", silent)
        
    if modified_distances:
        process.log("Using custom distances provided by the user.", silent)
        max_value = max(y for x in distances.values() for y in x)
        if max_value > 6:
            context.epsilon = max_value - 6

    process.log("\n", silent)
    process_func = process.single if core is None else process.multi_batch
    process_func(file_list, context)
    
    process.log("\n------------------------------------\n", silent)
    process.log(f"Total time elapsed: {(timer() - global_time_start):.3f}s\n", silent)


if __name__ == "__main__":
    main()
