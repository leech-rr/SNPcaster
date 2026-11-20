#!/usr/bin/env python3
"""
Comprehensive test data generator for snpeff_pipeline E2E testing.
Creates:
- test_data/ref.fasta (seqid: BA000007.3, 15350 bp)
- test_data/ref.gbk (LOCUS: BA000007.3, 14 CDSs, 15350 bp)
- test_data/ref.gff (GFF3 format for the 14 CDSs)
- test_data/ref.gtf (GTF format for the 14 CDSs)
- test_data/input_vcf/ (CHROM: BA000007.3, matching ref - for baseline tests)
- test_data/input_vcf_auto_map/ (CHROM: Sakai_BA000007, differing - for auto-map test)
- test_data/input_vcf_explicit_map/ (CHROM: chr_explicit, differing - for explicit map test)
- test_data/strain_list_*.txt
- test_data/chr_map_explicit.txt
- test_data/expected/
"""
import os
import re
import shutil
import subprocess
import sys

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_DIR = os.path.join(BASE_DIR, "ref-O157")
TEST_DATA_DIR = os.path.join(BASE_DIR, "test_data")
PIPELINE_SCRIPT = os.path.join(BASE_DIR, "snpeff_pipeline.py")

TARGET_LENGTH = 15350
CHR_REF = "BA000007.3"
CHR_AUTO = "Sakai_BA000007"
CHR_EXPLICIT = "chr_explicit"


def create_all():
    os.makedirs(TEST_DATA_DIR, exist_ok=True)

    # 1. Extract sequence
    src_fasta = os.path.join(REF_DIR, "Sakai-chromosome.fasta")
    seq_chunks = []
    total_len = 0
    with open(src_fasta, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            cleaned = line.strip()
            seq_chunks.append(cleaned)
            total_len += len(cleaned)
            if total_len >= TARGET_LENGTH:
                break
    full_seq = "".join(seq_chunks)[:TARGET_LENGTH]
    print(f"Extracted {len(full_seq)} bp sequence")

    # 2. Write ref.fasta
    ref_fasta = os.path.join(TEST_DATA_DIR, "ref.fasta")
    with open(ref_fasta, "w") as f:
        f.write(f">{CHR_REF}\n")
        for i in range(0, len(full_seq), 60):
            f.write(full_seq[i:i+60] + "\n")
    print(f"Created {ref_fasta}")

    # 3. Write ref.gbk
    src_gbk = os.path.join(REF_DIR, "Sakai-chromosome.gbk")
    ref_gbk = os.path.join(TEST_DATA_DIR, "ref.gbk")
    header_features = []
    with open(src_gbk, "r") as f:
        for line_num, line in enumerate(f, start=1):
            if line_num <= 313:
                if line.startswith("LOCUS"):
                    line = f"LOCUS       {CHR_REF:<22} {TARGET_LENGTH} bp    DNA     circular BCT 07-JUN-2018\n"
                elif "     source          1.." in line:
                    line = f"     source          1..{TARGET_LENGTH}\n"
                header_features.append(line)
            else:
                break

    with open(ref_gbk, "w") as f:
        for line in header_features:
            f.write(line)
        f.write("ORIGIN\n")
        for i in range(0, len(full_seq), 60):
            chunk = full_seq[i:i+60].lower()
            blocks = [chunk[j:j+10] for j in range(0, len(chunk), 10)]
            f.write(f"{i+1:9d} " + " ".join(blocks) + "\n")
        f.write("//\n")
    print(f"Created {ref_gbk}")

    # 4. Generate ref.gff using bp_genbank2gff3 and fix CDS phase via AGAT
    raw_gff = os.path.join(TEST_DATA_DIR, "ref.raw.gff")
    cmd = [
        "conda", "run", "--no-capture-output", "-n", "snippy",
        "perl", "-S", "bp_genbank2gff3",
        "-noinfer", "--filter", "gene", "--filter", "exon",
        "-out", "stdout", ref_gbk
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    gff_lines = []
    for line in res.stdout.splitlines(keepends=True):
        if line.startswith("##FASTA") or (gff_lines and gff_lines[-1].startswith("##FASTA")):
            break
        if line.startswith("##sequence-region"):
            line = f"##sequence-region {CHR_REF} 1 {TARGET_LENGTH}\n"
        elif line.startswith("# working on region:"):
            line = f"# working on region:{CHR_REF}\n"
        gff_lines.append(line)
    with open(raw_gff, "w") as f:
        f.writelines(gff_lines)

    ref_gff = os.path.join(TEST_DATA_DIR, "ref.gff")
    if os.path.exists(ref_gff):
        os.remove(ref_gff)
    agat_fix_cmd = [
        "mamba", "run", "--no-capture-output", "-n", "agat",
        "agat_sp_fix_cds_phases.pl",
        "--gff", raw_gff,
        "--fasta", ref_fasta,
        "-o", ref_gff
    ]
    subprocess.run(agat_fix_cmd, check=True, capture_output=True)
    if os.path.exists(raw_gff):
        os.remove(raw_gff)

    # Post-process ref.gff so that IDs match locus_tag
    with open(ref_gff, "r") as f:
        gff_lines = f.readlines()
    new_gff_lines = []
    for line in gff_lines:
        if line.startswith("#") or "\tregion\t" in line:
            new_gff_lines.append(line)
            continue
        parts = line.strip().split("\t")
        attrs_list = [x for x in parts[8].split(";") if "=" in x]
        attr_dict = dict(x.split("=", 1) for x in attrs_list)
        lt = attr_dict.get("locus_tag")
        ftype = parts[2]
        if lt:
            if ftype == "gene":
                attr_dict["ID"] = lt
            elif ftype == "mRNA":
                attr_dict["ID"] = lt
                attr_dict["Parent"] = lt
            elif ftype == "exon":
                attr_dict["ID"] = f"exon-{lt}"
                attr_dict["Parent"] = lt
            elif ftype == "CDS":
                attr_dict["ID"] = lt
                attr_dict["Parent"] = lt
        ordered_keys = [k for k in ["ID", "Parent", "Name"] if k in attr_dict]
        ordered_keys += [k for k in attr_dict if k not in ordered_keys]
        parts[8] = ";".join(f"{k}={attr_dict[k]}" for k in ordered_keys)
        new_gff_lines.append("\t".join(parts) + "\n")
    with open(ref_gff, "w") as f:
        f.writelines(new_gff_lines)
    print(f"Created {ref_gff} with fixed CDS phases and locus_tag IDs")

    # 5. Generate ref.gtf using AGAT and post-process IDs
    ref_gtf = os.path.join(TEST_DATA_DIR, "ref.gtf")
    if os.path.exists(ref_gtf):
        os.remove(ref_gtf)
    agat_gtf_cmd = [
        "mamba", "run", "--no-capture-output", "-n", "agat",
        "agat_convert_sp_gff2gtf.pl",
        "--gff", ref_gff,
        "-o", ref_gtf
    ]
    subprocess.run(agat_gtf_cmd, check=True, capture_output=True)

    with open(ref_gtf, "r") as f:
        gtf_lines = f.readlines()
    new_gtf_lines = []
    for line in gtf_lines:
        if line.startswith("#") or "\tregion\t" in line:
            new_gtf_lines.append(line)
            continue
        parts = line.strip().split("\t")
        attrs_str = parts[8].strip()
        m = re.search(r'locus_tag "([^"]+)"', attrs_str)
        if m:
            lt = m.group(1)
            attrs_str = re.sub(r'gene_id "[^"]+"', f'gene_id "{lt}"', attrs_str)
            attrs_str = re.sub(r'transcript_id "[^"]+"', f'transcript_id "{lt}"', attrs_str)
            attrs_str = re.sub(r'ID "[^"]+"', f'ID "{lt}"', attrs_str)
            attrs_str = re.sub(r'Parent "[^"]+"', f'Parent "{lt}"', attrs_str)
        parts[8] = attrs_str + (" " if not attrs_str.endswith(" ") else "")
        new_gtf_lines.append("\t".join(parts) + "\n")
    with open(ref_gtf, "w") as f:
        f.writelines(new_gtf_lines)
    print(f"Created {ref_gtf} using AGAT with locus_tag IDs")

    # 5.1 Generate ref_diff_locus_version.gbk (LOCUS != VERSION for testing)
    ref_diff = os.path.join(TEST_DATA_DIR, "ref_diff_locus_version.gbk")
    with open(ref_gbk, "r") as fin, open(ref_diff, "w") as fout:
        for line in fin:
            if line.startswith("LOCUS"):
                line = f"LOCUS       BA000007               {TARGET_LENGTH} bp    DNA     circular BCT 07-JUN-2018\n"
            fout.write(line)
    print(f"Created {ref_diff}")

    # 6. Create standard VCFs (CHROM: BA000007.3 - matches ref for standard pipeline tests)
    std_vcf_a = os.path.join(TEST_DATA_DIR, "input_vcf", "strain_A", "snps.vcf")
    std_vcf_b = os.path.join(TEST_DATA_DIR, "input_vcf", "strain_B", "snps.vcf")
    os.makedirs(os.path.dirname(std_vcf_a), exist_ok=True)
    os.makedirs(os.path.dirname(std_vcf_b), exist_ok=True)
    _filter_vcf(
        os.path.join(REF_DIR, "input-vcf", "DRR147067", "snps.vcf"),
        std_vcf_a,
        TARGET_LENGTH,
        rename_chr=CHR_REF
    )
    _filter_vcf(
        os.path.join(REF_DIR, "input-vcf", "SRR9694420", "snps.vcf"),
        std_vcf_b,
        TARGET_LENGTH,
        rename_chr=CHR_REF
    )
    print("Created standard input_vcf (matching CHROM)")

    # 7. Create differing VCF for Test 5: Auto-mapping (CHROM: Sakai_BA000007)
    auto_vcf_a = os.path.join(TEST_DATA_DIR, "input_vcf_auto_map", "strain_A", "snps.vcf")
    os.makedirs(os.path.dirname(auto_vcf_a), exist_ok=True)
    _filter_vcf(
        os.path.join(REF_DIR, "input-vcf", "DRR147067", "snps.vcf"),
        auto_vcf_a,
        TARGET_LENGTH,
        rename_chr=CHR_AUTO
    )
    print("Created input_vcf_auto_map (differing CHROM for auto-map test)")

    # 8. Create differing VCF for Test 6: Explicit --chr_map (CHROM: chr_explicit)
    explicit_vcf_a = os.path.join(TEST_DATA_DIR, "input_vcf_explicit_map", "strain_A", "snps.vcf")
    os.makedirs(os.path.dirname(explicit_vcf_a), exist_ok=True)
    _filter_vcf(
        os.path.join(REF_DIR, "input-vcf", "DRR147067", "snps.vcf"),
        explicit_vcf_a,
        TARGET_LENGTH,
        rename_chr=CHR_EXPLICIT
    )
    chr_map_explicit = os.path.join(TEST_DATA_DIR, "chr_map_explicit.txt")
    with open(chr_map_explicit, "w") as f:
        f.write(f"{CHR_REF}\t{CHR_EXPLICIT}\n")
    print("Created input_vcf_explicit_map and chr_map_explicit.txt")


    # 9. Create strain lists
    strain_single = os.path.join(TEST_DATA_DIR, "strain_list_single.txt")
    with open(strain_single, "w") as f:
        f.write("strain_A\n")
    with open(os.path.join(TEST_DATA_DIR, "strain_list_multi.txt"), "w") as f:
        f.write("strain_A\nstrain_B\n")
    with open(os.path.join(TEST_DATA_DIR, "strain_list_missing.txt"), "w") as f:
        f.write("strain_A\nnonexistent_strain\n")

    # 10. Generate expected folder using standard input_vcf
    print("Generating expected...")
    _generate_expected(ref_fasta, ref_gbk, "expected")
    # Clean up old separated dirs if they exist
    for old_dir in ["expected_gbk", "expected_gff", "expected_gtf"]:
        p = os.path.join(TEST_DATA_DIR, old_dir)
        if os.path.exists(p):
            shutil.rmtree(p, ignore_errors=True)
    print("All expected outputs generated successfully.")


def _filter_vcf(src_vcf: str, dst_vcf: str, max_pos: int, rename_chr: str = None):
    with open(src_vcf, "r") as fin, open(dst_vcf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                if rename_chr and line.startswith("##contig=<ID="):
                    line = re.sub(r"ID=[^,>]+", f"ID={rename_chr}", line)
                fout.write(line)
            else:
                parts = line.split("\t")
                pos = int(parts[1])
                if pos <= max_pos:
                    if rename_chr:
                        parts[0] = rename_chr
                        line = "\t".join(parts)
                    fout.write(line)
                else:
                    break


def _generate_expected(ref_fa: str, ann_file: str, expected_dirname: str):
    exp_dir = os.path.join(TEST_DATA_DIR, expected_dirname)
    db_tmp = os.path.join(TEST_DATA_DIR, f"tmp_db_{expected_dirname}")
    shutil.rmtree(exp_dir, ignore_errors=True)
    shutil.rmtree(db_tmp, ignore_errors=True)

    cmd = [
        sys.executable, PIPELINE_SCRIPT, "run",
        ref_fa, ann_file,
        os.path.join(TEST_DATA_DIR, "input_vcf"),
        os.path.join(TEST_DATA_DIR, "strain_list_multi.txt"),
        exp_dir,
        "--db_dir", db_tmp
    ]
    subprocess.run(cmd, check=True, capture_output=True)
    shutil.rmtree(db_tmp, ignore_errors=True)


if __name__ == "__main__":
    create_all()
