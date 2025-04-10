import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Create a summray for quality check result.')
parser.add_argument('coverage_result',      type=str,   help='Coverage result')
parser.add_argument('--checkm-result',      type=str,   help='Checkm result')
parser.add_argument('--assembly-summary',   type=str,   help='Assembly summary',)
parser.add_argument('--coverage-thresh',    type=float, help='Threshold for coverage')
parser.add_argument('--max-contigs',        type=int,   help='Maximum number of contigs')
parser.add_argument('--min-gen-size',       type=float, help='Minimum genome size')
parser.add_argument('--max-gen-size',       type=float, help='Maximum genome size')
parser.add_argument('--min-completeness',   type=float, help='Minimum completeness')
parser.add_argument('--max-contamination',  type=float, help='Maximum contamination')

args = parser.parse_args()

input_file_1 = args.coverage_result
input_file_2 = args.checkm_result
input_file_3 = args.assembly_summary
output_file = 'qc_results.xlsx'

cr_coverage=args.coverage_thresh
cr_max_contigs=args.max_contigs
cr_min_gen_size=args.min_gen_size
cr_max_gen_size=args.max_gen_size
cr_min_completeness=args.min_completeness
cr_max_contamination=args.max_contamination

print(f"Coverage file: {input_file_1}")
print(f"Checkm file: {input_file_2}")
print(f"Assembly summary file: {input_file_3}")
print(f"Output file: {output_file}")
print(f"Coverage threshold: {cr_coverage}")
print(f"Max contigs: {cr_max_contigs}")
print(f"Min genome size: {cr_min_gen_size}")
print(f"Max genome size: {cr_max_gen_size}")
print(f"Min completeness: {cr_min_completeness}")
print(f"Max contamination: {cr_max_contamination}")

df_coverage = pd.read_table(input_file_1)
if args.checkm_result is None:
    df_checkm = pd.DataFrame()
else:
    df_checkm = pd.read_table(input_file_2)
if args.assembly_summary is None:
    df_assembly = pd.DataFrame()
else:
    df_assembly = pd.read_table(input_file_3)

n=len(df_coverage)

cols = ['strain','no.of_reads','total_read_length','coverage','#contigs','largest_contig','total_length','marker_lineage','completeness','contamination', \
        'qc_results','qc_coverage','qc_#_contigs','qc_total_length','qc_completeness','qc_contamination']
df_qc = pd.DataFrame(index=[], columns=cols)

NA = 'N/A'
PASS = 'pass'
FAIL = 'fail'
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
        if(coverage >= float(cr_coverage)):
            qc_coverage=PASS
            qc_results = PASS if qc_results != FAIL else qc_results
        else:
            qc_coverage=FAIL
            qc_results = FAIL

    if(cr_max_contigs is not None):
        if(int(contigs) <= int(cr_max_contigs)):
            qc_contigs=PASS
            qc_results = PASS if qc_results != FAIL else qc_results
        else:
            qc_contigs=FAIL
            qc_results = FAIL

    if(cr_min_gen_size is not None \
       or cr_max_gen_size is not None):
        if(float(Total_length) >= float(cr_min_gen_size or 0.0)*10**6 \
           and float(Total_length) <= float(cr_max_gen_size or 10**6)*10**6):
            qc_Total_length=PASS
            qc_results = PASS if qc_results != FAIL else qc_results
        else:
            qc_Total_length=FAIL
            qc_results = FAIL

    if(cr_min_completeness is not None):
        if(float(Completeness) >= float(cr_min_completeness)):
            qc_completeness=PASS
            qc_results = PASS if qc_results != FAIL else qc_results
        else:
            qc_completeness=FAIL
            qc_results = FAIL

    if(cr_max_contamination is not None):
        if(float(Contamination) <= float(cr_max_contamination)):
            qc_contamination=PASS
            qc_results = PASS if qc_results != FAIL else qc_results
        else:
            qc_contamination=FAIL
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
    df_coverage.to_excel(writer, sheet_name="coverage", index=False)
    df_checkm.to_excel(writer, sheet_name="checkm_results", index=False)
    df_assembly.to_excel(writer, sheet_name="assembly_summary", index=False)
    df_qc.to_excel(writer,sheet_name="qc", index=False)

