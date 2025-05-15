"""
Author: Rafael Lemos - rafaellemos42@gmail.com
Date: 15/05/2025

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
            process_result(result, context.output)
        except Exception as e:
            print(f"Error: {e}")


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
        from psutil import Process, cpu_count
        
        if isinstance(core, list):
            Process(os.getpid()).cpu_affinity(core)
            if core[-1] - core[0] == len(core) - 1:
                print(f"Running on cores {core[0]} to {core[-1]}\nTotal number of cores: {len(core)}")
            else:
                print(f"Running on cores: {', '.join(map(str, core))}\nTotal number of cores: {len(core)}")
            num_cores = len(core)
        else:
            print(f"Running on {core} cores (automatically selected by OS)")
            num_cores = core
         
        batch_size = max(1, len(file_list) // num_cores)
        print(f"Number of files: {len(file_list)} | Batch size: {batch_size} files per core")
        print()
        
        with ProcessPoolExecutor(max_workers=num_cores) as executor:
            futures = {executor.submit(process_batch, batch, context): batch
                    for batch in batch_generator(file_list, batch_size)}

            for future in as_completed(futures):  # Wait for all batches to complete
                try:
                    future.result()  # Process results from batch
                except Exception as e:
                    print(f"Error processing batch: {e}")
                finally:
                    del futures[future] 
    except ImportError:
        print("Error.")
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
            process_result(result, context.output)


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
        parsed_data = parser.parse_pdb(file_path) if file_path.endswith(".pdb") else parser.parse_cif(file_path)

        if parsed_data.true_count() > 10000:  # Skip very large proteins (customizable)
            print(f"Skipping ID '{parsed_data.id}'. Size: {parsed_data.true_count()} residues") 
            if context.output:
                with open(f"{context.output}/big.csv", "a") as f:
                    f.write(f"{parsed_data.id},{parsed_data.title},{parsed_data.true_count()},x\n")
            return None

        contacts_list = contacts.contact_detection(parsed_data, context.custom_distances, context.epsilon)
        process_time = timer() - start_time
        return parsed_data, contacts_list, process_time

    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return None


def process_result(result, output):
    """
    Handles the result of processing a file.

    Args:
        result (tuple): A tuple containing the processed Protein object, contacts list, and processing time.
        output (str): The directory where output files will be saved.
    """
    if result:
        protein, contacts_list, process_time = result
        output_data = f"ID: {protein.id} | Size: {protein.true_count():<7} | Contacts: {len(contacts_list):<7} | Time: {process_time:.3f}s"
        print(output_data)
        
        if output:
            output_folder = f"{output}/{protein.id}/"
            
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
            
            with open(f"{output_folder}/{protein.id}_contacts.csv","w") as f:
                f.write(contacts.show_contacts(contacts_list))


def validate_categories(categories):
    for key, (min_val, max_val) in categories.items():
        if min_val < 0 or max_val < 0:
            raise ValueError(f"Invalid values for '{key}': values must be positive.")
        if min_val >= max_val:
            raise ValueError(f"Invalid range for '{key}': min ({min_val}) must be less than max ({max_val}).")
    return categories