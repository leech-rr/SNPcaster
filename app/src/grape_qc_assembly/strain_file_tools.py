from file_tools import read_tsv_file

def extract_strains(fastq_list_path: str) -> tuple[str, ...]:
    lines = read_tsv_file(fastq_list_path)
    return tuple(
        line[0] for line in lines
    )

