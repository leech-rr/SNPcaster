import pandas as pd
import argparse
import os
import subprocess
import logging

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description='Create a summary for quality check result.')
parser.add_argument('--checkm-result',      type=str,   help='Checkm result file', required=True)
parser.add_argument('--assembly-summary',   type=str,   help='Assembly summary file', required=True)
parser.add_argument('--assembly-list',      type=str,   help='Assembly list file mapping strains to fastq files', required=True)
parser.add_argument('--fastq-dir',         type=str,   help='Directory containing fastq files for coverage calculation', required=True)
parser.add_argument('--threads',           type=int,   default=1, help='Number of threads for seqkit')
parser.add_argument('--coverage-thresh',    type=float, help='Threshold for coverage')
parser.add_argument('--max-contigs',        type=int,   help='Maximum number of contigs')
parser.add_argument('--min-gen-size',       type=float, help='Minimum genome size')
parser.add_argument('--max-gen-size',       type=float, help='Maximum genome size')
parser.add_argument('--min-completeness',   type=float, help='Minimum completeness')
parser.add_argument('--max-contamination',  type=float, help='Maximum contamination')

args = parser.parse_args()

def calculate_read_stats(fastq_dir, fastq_filenames, threads=1):
    """Calculate total read length for a strain using seqkit stats and specific files."""
    fastq_files = []
    for fname in fastq_filenames:
        if fname.strip():
            fastq_files.append(os.path.join(fastq_dir, fname.strip()))
    
    if not fastq_files:
        logger.warning(f"No fastq files provided to calculate_read_stats in {fastq_dir}")
        return 0, 0

    # Ensure files exist before passing to seqkit
    valid_files = [f for f in fastq_files if os.path.exists(f)]
    if not valid_files:
        logger.warning(f"None of the provided fastq files exist in {fastq_dir}: {fastq_files}")
        return 0, 0

    total_reads = 0
    total_len = 0
    
    try:
        # Use conda run -n seqkit seqkit stats if possible, or just seqkit if in environment
        # Based on previous implementation, we try conda run
        cmd = ['conda', 'run', '-n', 'seqkit', 'seqkit', 'stats', '-T', '-j', str(threads)] + valid_files
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        lines = result.stdout.strip().split('\n')
        if len(lines) < 2:
            return 0,0
            
        header = lines[0].split('\t')
        num_seqs_idx = header.index('num_seqs')
        sum_len_idx = header.index('sum_len')
        
        for line in lines[1:]:
            parts = line.split('\t')
            total_reads += int(parts[num_seqs_idx].replace(',', ''))
            total_len += int(parts[sum_len_idx].replace(',', ''))
            
    except Exception as e:
        logger.error(f"Error running seqkit for files {valid_files}: {e}")
        return 0, 0
        
    return total_reads, total_len

NA = 'N/A'
PASS = 'pass'
FAIL = 'fail'

def to_float_or_none(x):
    """Return float(x) or None if conversion fails or x is NA-like."""
    try:
        if x is None:
            return None
        sx = str(x).strip()
        if sx == '' or sx.upper() == NA:
            return None
        # remove commas commonly used in thousands separators
        sx = sx.replace(',', '')
        return float(sx)
    except Exception as e:
        # Raise a clear ValueError so caller knows which value failed to parse
        raise ValueError(f"to_float_or_none: cannot convert {x!r} to float") from e

def to_int_or_none(x):
    """Return int(x) or None if conversion fails or x is NA-like."""
    fv = to_float_or_none(x)
    if fv is None:
        return None
    try:
        return int(fv)
    except Exception as e:
        # Raise a clear ValueError so caller knows which value failed to parse
        raise ValueError(f"to_int_or_none: cannot convert {x!r} to int") from e

input_file_2 = args.checkm_result
input_file_3 = args.assembly_summary
output_file = 'qc_results.xlsx'

cr_coverage=args.coverage_thresh
cr_max_contigs=args.max_contigs
cr_min_gen_size=args.min_gen_size
cr_max_gen_size=args.max_gen_size
cr_min_completeness=args.min_completeness
cr_max_contamination=args.max_contamination

print(f"Checkm file: {input_file_2}")
print(f"Assembly summary file: {input_file_3}")
print(f"Assembly list file: {args.assembly_list}")
print(f"Fastq directory: {args.fastq_dir}")
print(f"Output file: {output_file}")
print(f"Coverage threshold: {cr_coverage}")
print(f"Max contigs: {cr_max_contigs}")
print(f"Min genome size: {cr_min_gen_size}")
print(f"Max genome size: {cr_max_gen_size}")
print(f"Min completeness: {cr_min_completeness}")
print(f"Max contamination: {cr_max_contamination}")

df_assembly = pd.read_table(input_file_3)

# Parse assembly list mapping
strain_to_fastq = {}
with open(args.assembly_list, 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) > 1:
            strain_to_fastq[parts[0]] = parts[1:]
        else:
            strain_to_fastq[parts[0]] = []

# Build df_coverage from assembly summary and fastq stats
coverage_data = []
for i in range(len(df_assembly)):
    strain = str(df_assembly.iloc[i, 0]).strip()
    if strain.endswith('_s'): # skip scaffolds
        continue
        
    total_len_assembly = to_float_or_none(df_assembly.iloc[i, 15]) # 'Total length' column
    
    if strain not in strain_to_fastq:
        raise KeyError(f"Strain '{strain}' found in assembly summary but missing from assembly list.")
    reads, total_read_len = calculate_read_stats(args.fastq_dir, strain_to_fastq[strain], args.threads)
        
    coverage = total_read_len / total_len_assembly if total_len_assembly and total_len_assembly > 0 else 0
    coverage_data.append([strain, reads, total_read_len, coverage])

df_coverage = pd.DataFrame(coverage_data, columns=['strain', 'no. of reads', 'total read length', 'coverage'])
df_checkm = pd.read_table(input_file_2)

n=len(df_coverage)

cols = ['strain','no.of_reads','total_read_length','coverage','#contigs','largest_contig','total_length','marker_lineage','completeness','contamination', \
        'qc_results','qc_coverage','qc_#_contigs','qc_total_length','qc_completeness','qc_contamination']
df_qc = pd.DataFrame(index=[], columns=cols)

for i in range(len(df_coverage)):
    strain_cov = df_coverage.iloc[i,0]
    no_of_reads= df_coverage.iloc[i,1]
    total_read_length=df_coverage.iloc[i,2]
    coverage=df_coverage.iloc[i,3]

    contigs=NA
    Largest_contig=NA
    Total_length=NA 
    Marker_lineage=NA
    Completeness=NA
    Contamination=NA
    qc_coverage=NA
    qc_contigs=NA
    qc_Total_length=NA
    qc_completeness=NA
    qc_contamination=NA
    qc_results=NA
         
    for j in range(len(df_checkm)):
        if(strain_cov == df_checkm.iloc[j,0]):
            Marker_lineage=df_checkm.iloc[j,1]
            Completeness=df_checkm.iloc[j,11]
            Contamination=df_checkm.iloc[j,12]
    for k in range(len(df_assembly)):
        if(strain_cov == df_assembly.iloc[k,0]):
            contigs=df_assembly.iloc[k,13]
            Largest_contig=df_assembly.iloc[k,14]
            Total_length=df_assembly.iloc[k,15]           
 
    if(cr_coverage is not None):
        # coverage from df_coverage may be non-numeric; parse safely
        coverage_f = to_float_or_none(coverage)
        if coverage_f is not None and coverage_f >= cr_coverage:
            qc_coverage = PASS
            qc_results = PASS if qc_results != FAIL else qc_results
        else:
            qc_coverage = FAIL
            qc_results = FAIL

    if(cr_max_contigs is not None):
        contigs_i = to_int_or_none(contigs)
        if contigs_i is not None and contigs_i <= cr_max_contigs:
            qc_contigs = PASS
            qc_results = PASS if qc_results != FAIL else qc_results
        else:
            qc_contigs = FAIL
            qc_results = FAIL

    if(cr_min_gen_size is not None \
       or cr_max_gen_size is not None):
        total_len_f = to_float_or_none(Total_length)
        min_size = (cr_min_gen_size or 0.0)
        max_size = (cr_max_gen_size or 10**6)
        # convert genome sizes to bases (user supplies Mb)
        if total_len_f is not None and total_len_f >= (min_size * 10**6) \
           and total_len_f <= (max_size * 10**6):
            qc_Total_length = PASS
            qc_results = PASS if qc_results != FAIL else qc_results
        else:
            qc_Total_length = FAIL
            qc_results = FAIL

    if(cr_min_completeness is not None):
        completeness_f = to_float_or_none(Completeness)
        if completeness_f is not None and completeness_f >= cr_min_completeness:
            qc_completeness = PASS
            qc_results = PASS if qc_results != FAIL else qc_results
        else:
            qc_completeness = FAIL
            qc_results = FAIL

    if(cr_max_contamination is not None):
        contamination_f = to_float_or_none(Contamination)
        if contamination_f is not None and contamination_f <= cr_max_contamination:
            qc_contamination = PASS
            qc_results = PASS if qc_results != FAIL else qc_results
        else:
            qc_contamination = FAIL
            qc_results = FAIL

    # convert to upper case
    qc_results=str.upper(qc_results)

    data_new=[[strain_cov,no_of_reads,total_read_length,coverage,contigs,Largest_contig,Total_length,Marker_lineage,Completeness,Contamination,\
              qc_results, qc_coverage, qc_contigs,qc_Total_length,qc_completeness,qc_contamination]]
    df_new=pd.DataFrame(data_new,columns=cols)
    df_qc=pd.concat([df_qc,df_new],axis=0)

print(df_qc)

df_qc.to_csv('qc_results.tsv', sep='\t', index=False)

with pd.ExcelWriter(output_file) as writer:
    df_checkm.to_excel(writer, sheet_name="checkm_results", index=False)
    df_assembly.to_excel(writer, sheet_name="assembly_summary", index=False)
    df_qc.to_excel(writer,sheet_name="qc", index=False)

