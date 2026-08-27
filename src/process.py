"""
COCaDA v1.6
Author: Rafael Lemos - rafaellemos@ufmg.br
Date: 26/03/2026

License: MIT License
"""

import os
from itertools import islice
from timeit import default_timer as timer
import src.contacts as contacts
import src.parser as parser

def single(file_list, context):
    """
    Processes a list of files in single-core mode.

    Args:
        file_list (list): List of file paths to process.
        context (ProcessingContext): Context object containing parameters such as core, output, and region.

    This function processes each file in the list sequentially, detects contacts, and outputs the results to the console or to a file, depending on the 'output' flag.
    """
    for file in file_list:
        try:
            result = process_file(file, context)
            process_result(result, context)
        except Exception as e:
            log(f"Error: {e}")


def multi_batch(file_list, context):
    """
    Distributes the processing of files across multiple cores in batches.

    Args:
        file_list (list): List of file paths to process.
        context (ProcessingContext): Context object containing parameters such as core, output, and region.
    """
    core = context.core

    try:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        from psutil import Process
        
        if isinstance(core, list):
            try:
                Process(os.getpid()).cpu_affinity(core)
            except AttributeError:
                log("Setting core affinity is not supported on this OS. Using cores automatically selected by OS", context.silent)
                
            if core[-1] - core[0] == len(core) - 1:
                log(f"Running on cores {core[0]} to {core[-1]}\nTotal number of cores: {len(core)}", context.silent)
            else:
                log(f"Running on cores: {', '.join(map(str, core))}\nTotal number of cores: {len(core)}", context.silent)
            num_cores = len(core)
        else:
            log(f"Running on {core} cores (automatically selected by OS)", context.silent)
            num_cores = core
         
        batch_size = max(1, len(file_list) // num_cores)
        log(f"Number of files: {len(file_list)} | Batch size: {batch_size} files per core", context.silent)
        log("\n", context.silent)
        
        with ProcessPoolExecutor(max_workers=num_cores) as executor:
            futures = {executor.submit(process_batch, batch, context): batch
                    for batch in batch_generator(file_list, batch_size)}

            for future in as_completed(futures):  # Wait for all batches to complete
                try:
                    future.result()  # Process results from batch
                except Exception as e:
                    log(f"Error processing batch: {e}")
                finally:
                    del futures[future] 
    except ImportError:
        log("Error.")
        exit(1)


def process_batch(batch, context):
    """
    Processes a single batch of files sequentially.

    Args:
        batch (list): List of file paths in the batch.
        context (ProcessingContext): Context object containing parameters such as core, output, and region.
    """
    for file_path in batch:
        result = process_file(file_path, context)
        if result:
            process_result(result, context)


def batch_generator(file_list, batch_size):
    """
    Generates batches from the file list.

    Args:
        file_list (list): List of file paths to split into batches.
        batch_size (int): Maximum number of files per batch.

    Yields:
        list of str: Next batch of file paths.
    """
    it = iter(file_list)
    while batch := list(islice(it, batch_size)):
        yield batch


def process_file(file_path, context):
    """
    Processes a single file for contact detection.

    Args:
        file_path (str): Path to the file to be processed.
        context (ProcessingContext): Context object containing parameters such as core, output, and region.

    Returns:
        tuple: A tuple containing the processed Protein object, the list of detected contacts, and the processing time.
        None: If the file cannot be processed or an error occurs.

    This function parses the PDB or mmCIF file, detects contacts, and returns the results. If an error occurs during processing, it logs the error and returns None.
    """
    start_time = timer()

    try:
        parsed_data, ph = parser.parse_pdb(file_path) if file_path.endswith(".pdb") else parser.parse_cif(file_path, context.propedia)
                
        if context.ph is None:
            uncertainty_flags, local_contact_types = contacts.change_protonation(ph, context.silent)
            if ph != 7.4:
                log(f"Found experimental protein pH value at {ph}. You can change this using the -ph flag.", context.silent)
                log(f"Changing protonation states of pH-sensitive atoms using pH value of {ph}.", context.silent)
            else:
                log("Defaulting pH value to 7.4.", context.silent)
        else:
            if context.ph == -1:
                context.ph = ph
            uncertainty_flags, local_contact_types = contacts.change_protonation(context.ph, context.silent)
            
        if parsed_data.true_count() > 10000:  # Skip very large proteins (customizable)
            log(f"Skipping ID '{parsed_data.id}'. Size: {parsed_data.true_count()} residues.\n", context.silent) 
            if context.output:
                with open(f"{context.output}/big.csv", "a") as f:
                    f.write(f"{parsed_data.id},{parsed_data.title},{parsed_data.true_count()},x,{ph}\n")
            return None
            
        contacts_list, interface_res, count_types = contacts.contact_detection(parsed_data, context, uncertainty_flags, local_contact_types)
        process_time = timer() - start_time

        return parsed_data, contacts_list, process_time, interface_res, count_types, ph

    except Exception as e:
        log(f"Error processing {file_path}: {e}")
        return None


def process_result(result, context):
    """
    Handles the result of processing a file.

    Args:
        result (tuple): A tuple containing the processed Protein object, contacts list, and processing time.
        output (str): The directory where output files will be saved.
    """
    if result:
        protein, contacts_list, process_time, interface_res, count_types, ph = result
        output, silent, web = context.output, context.silent, context.web
        ph = ph if context.ph is None else context.ph
        
        output_data = f"ID: {protein.id} | Size: {protein.true_count():<7} | Contacts: {len(contacts_list):<7} | pH: {ph:.2f} | Time: {process_time:.3f}s"
        count = '; '.join(f"{v[0]}: {v[1]:>5}" for v in count_types.values())
        log(output_data)
        log(f"{count}\n", silent)
                
        if output:
            os.makedirs(output, exist_ok=True)

            if web == 0: # default
                contacts_file = f"{output}/{protein.id}_contacts.csv"
            elif web == 1: # web
                contacts_file = f"{output}/contacts.csv"
                info_file = f"{output}/info.csv"
            elif web == 2: # database update
                base = f"{output}/{protein.id}"
                os.makedirs(base, exist_ok=True)
                contacts_file = f"{base}/{protein.id}_contacts.csv"
                info_file = f"{base}/{protein.id}_info.csv"
            
            with open(contacts_file,"w") as f:
                f.write(contacts.show_contacts(contacts_list))
                    
            # COCaDA-web exclusive
            if web == 1 or web == 2: # web or database update
                number_contacts = contacts.count_contacts(contacts_list)
                number_contacts = ','.join(map(str, number_contacts))
                with open(info_file,"w") as f:
                    f.write(f"{protein.id},{protein.title},{protein.true_count()},{len(contacts_list)},{number_contacts},{ph}")
            if web == 0 or web == 2: # default and database update
                list_file = f"{output}/list.csv"
                with open(list_file,"a") as f:
                    f.write(f"{protein.id},{protein.title},{protein.true_count()},{len(contacts_list)},{ph}\n")

def log(message, silent=False):
    if not silent:
        print(message)
