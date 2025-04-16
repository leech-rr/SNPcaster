#!/usr/bin/env python3

import os
import glob
import argparse
from posixpath import basename
import logging

from file_tools import read_file, write_tsv_file

logger = logging.getLogger()

def find_fastq_pairs(
    strain_list_file: str,
    target_folder: str,
    file_extension: str,
    paired_file: str,
    unpaired_file: str,
) -> None:
    """
    This function searches for paired fastq files for each strain in the given strain list,
    and outputs the results to a specified file.

    Args:
    strain_list_file (str): Path to the file containing the list of strain names.
    target_folder (str): Folder to search for files.
    file_extension (str): Extension of the files to search for.
    output_file (str): Path to the output file.
    """
    # Read strain names from the input file
    strain_list = read_file(strain_list_file)

    # Sort strain names by length to avoid partial matches (longest first)
    length_sorted_strain_list = sorted(strain_list, key=len, reverse=True)

    # Track used files to prevent reusing them
    used_files = set()

    # Lists to store paired and unpaired strains
    paired_strain_to_files = {}
    unpaired_strain_to_files = {}

    # Search for paired files
    for strain in length_sorted_strain_list:
        if strain in paired_strain_to_files \
            or strain in unpaired_strain_to_files:
            raise Exception(f"Duplicate strain name [{strain}] found in the list.")

        search_pattern = os.path.join(target_folder, f"{strain}*.{file_extension}")
        files = list(map(lambda p: basename(p), glob.glob(search_pattern)))

        # Filter out files that have already been used
        files = [f for f in files if f not in used_files]

        # If exactly two files are found, they are considered as a pair
        if len(files) == 2:
            paired_strain_to_files[strain] = sorted(files)
            used_files.update(files)
        else:
            # If more than one file is found but not a pair, list all the files
            unpaired_strain_to_files[strain] = sorted(files)

    # Convert to list of tuples for easier writing
    paired_strains = [
        (strain, *(paired_strain_to_files[strain])) 
        for strain in strain_list
        if strain in paired_strain_to_files
    ]
    unpaired_strains = [
        (strain, *unpaired_strain_to_files[strain])
        for strain in strain_list
        if strain in unpaired_strain_to_files
    ]

    # Write results to the output file
    write_tsv_file(paired_file, paired_strains)
    if len(unpaired_strains) > 0:
        write_tsv_file(unpaired_file, unpaired_strains)
        logger.warning(f"Unpaired files found for {len(unpaired_strains)} strain(s). Check the output file [{unpaired_file}].")

def main() -> None:
    """
    Main function to handle command-line arguments and call the find_fastq_pairs function.
    """

    try:
        # Set up command-line argument parsing
        parser = argparse.ArgumentParser(description="Find paired fastq files for strains listed in a file.")
        parser.add_argument('strain_list_file', type=str, help="Path to the file containing the list of strain names.")
        parser.add_argument('--target_folder', type=str, default=".", help="Folder to search for files (default: current directory).")
        parser.add_argument('--file_extension', type=str, default="fastq.gz", help="File extension to search for (default: fastq.gz).")
        parser.add_argument('--paired_list', type=str, default="list_fastq.tsv", help="Path to the file for strain paired files (default: list_fastq.tsv).")
        parser.add_argument('--unpaired_list', type=str, default="unpaired_fastq.tsv", help="Path to the file for strain unpaired files (default: unpaired_fastq.tsv).")

        # Parse command-line arguments
        args = parser.parse_args()

        # Call the find_fastq_pairs function with the parsed arguments
        find_fastq_pairs(
            strain_list_file=args.strain_list_file,
            target_folder=args.target_folder,
            file_extension=args.file_extension,
            paired_file=args.paired_list,
            unpaired_file=args.unpaired_list,
        )
    except Exception:
        logger.exception(f"Exit with following error...")
        exit(1)

if __name__ == "__main__":
    main()
