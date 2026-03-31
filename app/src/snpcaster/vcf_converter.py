#!/usr/bin/env python3
import argparse
import subprocess
import tempfile
import os
import shutil
import sys
import datetime

def parse_args():
    parser = argparse.ArgumentParser(description="Convert BactSNP output to well-formed VCF using bcftools.")
    parser.add_argument("source_fasta", help="BactSNP alignment FASTA")
    parser.add_argument("ref_fasta", help="Reference genome FASTA")
    parser.add_argument("target_id", help="Target sample ID to extract")
    return parser.parse_args()

def get_snpcaster_version():
    res = subprocess.run(["snpcaster.sh", "-v"], capture_output=True, text=True)
    if res.returncode == 0:
        for line in res.stdout.strip().split('\n'):
            if line.startswith("SNPCaster"):
                return escape_for_vcf(line)
    raise Exception("Could not determine SNPcaster version because `snpcaster.sh -v` failed.")

def escape_for_vcf(text):
    return text.replace(" ", "_").replace("\"", "\\\"") 

def run_command(cmd, check=True, use_snippy=True):
    full_cmd = ["conda", "run", "-n", "snippy"] + cmd if use_snippy else cmd
    print(f"Running: {' '.join(full_cmd)}", file=sys.stderr)
    result = subprocess.run(full_cmd, capture_output=True, text=True)
    if check and result.returncode != 0:
        print(f"Error running command: {' '.join(full_cmd)}\n{result.stderr}", file=sys.stderr)
        sys.exit(result.returncode)
    return result

def get_fasta_contig_info(fasta_path):
    contig_name = None
    contig_length = 0
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if not contig_name:
                    contig_name = line.strip()[1:].split()[0]
            elif contig_name:
                contig_length += len(line.strip())
    return contig_name, contig_length

def main():
    args = parse_args()

    contig_name, contig_length = get_fasta_contig_info(args.ref_fasta)
    if not contig_name:
        print("Failed to read contig name from reference FASTA", file=sys.stderr)
        sys.exit(1)

    # Create a temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_combined_fasta = os.path.join(temp_dir, "combined.fasta")
        
        # Step 1: Extraction
        try:
            with open(temp_combined_fasta, "w") as out_f:
                with open(args.ref_fasta, "r") as ref_f:
                    shutil.copyfileobj(ref_f, out_f)
                
                with open(args.source_fasta, "r") as src_f:
                    for line in src_f:
                        if line.startswith(">"):
                            out_f.write(f">{args.target_id}\n")
                        else:
                            out_f.write(line)
        except Exception as e:
            print(f"Failed to create combined FASTA: {e}", file=sys.stderr)
            sys.exit(1)

        raw_vcf = os.path.join(temp_dir, "raw.vcf")
        run_command(["snp-sites", "-v", temp_combined_fasta, "-o", raw_vcf])

        # Step 2: Reformatting and Reheadering via bcftools        

        # Generate custom headers to inject via bcftools
        file_date = datetime.date.today().strftime("%Y%m%d")
        custom_headers = os.path.join(temp_dir, "custom_headers.txt")
        with open(custom_headers, "w") as f_out:
            f_out.write(f'##fileDate={file_date}\n')
            f_out.write(f'##source={get_snpcaster_version()}\n')
            f_out.write(f'##reference={escape_for_vcf(os.path.basename(args.ref_fasta))}\n')
            
        # We need to rename the chromosome '1' (from snp-sites) to the actual contig name
        chr_rename_txt = os.path.join(temp_dir, "chr_rename.txt")
        with open(chr_rename_txt, "w") as f:
            f.write(f"1 {contig_name}\n")
            
        renamed_vcf = os.path.join(temp_dir, "renamed.vcf")
        run_command(["bcftools", "annotate", "--no-version", "-h", custom_headers, "--rename-chrs", chr_rename_txt, raw_vcf, "-o", renamed_vcf])
        
        sample_vcf = os.path.join(temp_dir, "sample.vcf")
        run_command(["bcftools", "view", "--no-version", "-s", args.target_id, renamed_vcf, "-o", sample_vcf])

        # Step 3: Ploidy and Normalization
        # bcftools +fixploidy sample.vcf -- -p ploidy_map.txt
        # If the plugin is not installed/working, it might fail. Let's try it.
        # Format for ploidy_map: CHROM FROM TO SEX PLOIDY
        # We assume 1 for haploid or just 2 for diploid if we want to force diploid
        ploidy_txt = os.path.join(temp_dir, "ploidy.txt")
        with open(ploidy_txt, "w") as f:
            f.write(f"{contig_name}\t1\t-\tM\t2\n")
            
        diploid_vcf = os.path.join(temp_dir, "diploid.vcf")
        fixploidy_result = run_command(["bcftools", "+fixploidy", "--no-version", sample_vcf, "-o", diploid_vcf, "--", "-p", ploidy_txt], check=False)
        
        if fixploidy_result.returncode != 0:
            print(f"Note: bcftools +fixploidy failed or plugin not found: {fixploidy_result.stderr}", file=sys.stderr)
            print("Falling back to manual text manipulation for haploid -> diploid conversion.", file=sys.stderr)
            
            # Simple fallback for GT 0 -> 0/0 and 1 -> 1/1
            with open(sample_vcf, 'r') as infile, open(diploid_vcf, 'w') as outfile:
                for line in infile:
                    if line.startswith('#'):
                        outfile.write(line)
                    else:
                        parts = line.strip().split('\t')
                        if len(parts) >= 10:
                            gt_info = parts[9].split(':')
                            if gt_info[0] == '0':
                                gt_info[0] = '0/0'
                            elif gt_info[0] == '1':
                                gt_info[0] = '1/1'
                            elif gt_info[0] == '.':
                                gt_info[0] = './.'
                            parts[9] = ':'.join(gt_info)
                        outfile.write('\t'.join(parts) + '\n')

        # Run bcftools norm
        final_vcf = f"{args.target_id}.vcf"
        run_command(["bcftools", "norm", "--no-version", "-f", args.ref_fasta, "-m", "-any", diploid_vcf, "-o", final_vcf])

        # bcftools keeps VCFv4.1 from snp-sites. Force it to VCFv4.2 via command line string replacement.
        subprocess.run(["sed", "-i", "s/^##fileformat=VCFv4.1/##fileformat=VCFv4.2/", final_vcf])

        print(f"Successfully created {final_vcf}", file=sys.stderr)

if __name__ == "__main__":
    main()
