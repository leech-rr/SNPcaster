import argparse
from argparse import Namespace
import logging
import os

from file_tools import read_file, write_file
from strain_file_tools import extract_strains

logger = logging.getLogger()


def parse_args() -> Namespace:
    parser = argparse.ArgumentParser(description='Create a summray of quast results for assemby.')
    parser.add_argument('quast_folder',      type=str, help='Input quast result folder')
    parser.add_argument('assembly_list',     type=str, help='Input assembly list')
    parser.add_argument('assembly_version',  type=str, help='Assembly version text')
    parser.add_argument('summary_file_path', type=str, help='Output summary file path')
    args = parser.parse_args()
    return args


def extract_quast_result(quast_result_file: str) -> tuple[str, str]:
    headers = []
    values = []
    quast_result = read_file(quast_result_file)
    for res in quast_result:
        cols = res.split("\t")
        headers.append(cols[0])
        values.append(cols[1])
    return "\t".join(headers), "\t".join(values)

def main():
    try:
        args = parse_args()
        assembly_list = args.assembly_list
        assembly_version = args.assembly_version
        summary_file_path = args.summary_file_path
        quast_folder = args.quast_folder
        strains = extract_strains(assembly_list)

        header_added = False
        merged_quast_result = []
        for strain in strains:
            quast_result_file = f"{quast_folder}/{strain}/report.tsv"
            header, result = extract_quast_result(quast_result_file)
            if not header_added:
                merged_quast_result.append(f"{header}\tAssembler version")
                header_added = True
            merged_quast_result.append(f"{result}\t{assembly_version}")

            # Scaffold result
            quast_result_file_scaffold = f"{quast_folder}/{strain}_s/report.tsv"
            if not os.path.exists(quast_result_file_scaffold):
                continue
            header, result = extract_quast_result(quast_result_file_scaffold)
            merged_quast_result.append(f"{result}\t{assembly_version}")
        write_file(summary_file_path, merged_quast_result)
    except Exception as e:
        logger.exception(f"Exit with following error...")
        exit(1)

if __name__ == "__main__":
    main()

