#!/bin/bash

# Function to extract a specific column and output it to stdout
# Example usage:
# extract_column_to_stdout "input.tsv" 1
extract_column_to_stdout() {
    input_file="$1"    # Input TSV file path
    column_number="$2" # Column number to extract

    if [[ ! -r "$input_file" ]]; then
        echo "Error: Input file '$input_file' can not read." >&2
        return 1
    fi
    # Check if the column number is specified and valid
    if [[ -z "$column_number" || "$column_number" -lt 1 ]]; then
        echo "Please specify a valid column number (an integer greater than or equal to 1)." >&2
        return 1
    fi
    awk -v col="$column_number" -F'\t' '{print $col}' "$input_file"
}

# Function to extract a specific column and output it to stdout OR a specified file path
# Example usage:
# extract_column "input.tsv" 1 "output.txt" 
extract_column() {
    input_file="$1"    # Input TSV file path
    column_number="$2" # Column number to extract
    output_file="${3:-}"   # Output file path (optional, if not specified, output to stdout)

    # Extract the specified column
    if [[ -z "$output_file" ]]; then
        extract_column_to_stdout "$input_file" "$column_number"        
    else
        extract_column_to_stdout "$input_file" "$column_number" > "$output_file"
    fi

    # Return the exit status of the extraction function
    return $?
}


# Function to extract strain names from fastq list TSV file
# Example usage:
# extract_column "list_fastq.tsv" "list.txt"
extract_strain_names() {
    input_file="$1"    # Input TSV file path
    output_file="${2:-}"   # (optional) Output file path (if not specified, output to stdout)
    extract_column "$input_file" 1 "$output_file"

    # Return the exit status of the extraction function
    return $?
}

# Function to extract strain names from fastq list TSV file
# Example usage:
# extract_column "list_fastq.tsv" "list.txt"
filter_fastq_list() {
    fastq_list="$1"    # Input fastq file path
    filter_list="$2"   # Strain list to filter

    # Check if the fastq list file exists
    if [[ ! -r "$fastq_list" ]]; then
        echo "Error: Fastq list file '$fastq_list' can not read." >&2
        return 1
    fi
    # Check if the filter list file exists
    if [[ ! -r "$filter_list" ]]; then
        echo "Error: Filter list file '$filter_list' can not read." >&2
        return 1
    fi
    # Extract the strain names from the filter list
    awk 'NR==FNR {stock[$1]; next} $1 in stock' "${filter_list}" "${fastq_list}"
}
