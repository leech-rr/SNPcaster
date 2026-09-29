#!/usr/bin/env python3
import argparse
import subprocess
import os
import sys
import shutil
import logging
import re
from dataclasses import dataclass, field
from typing import List, Optional, IO, Dict

# ==============================================================================
# Constants
# ==============================================================================
CUSTOM_DB_NAME = "snpeff-custom"
CONFIG_FILE_NAME = "snpEff.config"
SEQUENCES_FA_NAME = "sequences.fa"
GENES_GBK_NAME = "genes.gbk"
GENES_GFF_NAME = "genes.gff"
GENES_GTF_NAME = "genes.gtf"
DEFAULT_DB_DIR = "./custom_db"
DEFAULT_CODON_TABLE = "[Bacteria]"
DEFAULT_UD = 300

# Configure logger
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("snpeff_pipeline")


# ==============================================================================
# Parameter / Configuration Data Classes
# ==============================================================================
@dataclass
class EnvironmentConfig:
    """Configuration for Conda environments."""
    snpeff_env: str = "snpeff"
    bcftools_env: str = "snippy"
    snippy_env: str = "snippy"
    bioperl_env: str = "perl-bioperl-core"

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> "EnvironmentConfig":
        return cls(
            snpeff_env=args.snpeff_env,
            bcftools_env=args.bcftools_env,
            snippy_env=args.snippy_env,
            bioperl_env=getattr(args, "bioperl_env", "perl-bioperl-core"),
        )


@dataclass
class CustomDatabase:
    """Manages paths and parameters for a custom SnpEff database."""
    db_dir: str
    codon_table: str
    db_name: str = CUSTOM_DB_NAME
    config_filename: str = CONFIG_FILE_NAME
    sequences_fa_filename: str = SEQUENCES_FA_NAME
    genes_gbk_filename: str = GENES_GBK_NAME
    genes_gff_filename: str = GENES_GFF_NAME
    genes_gtf_filename: str = GENES_GTF_NAME

    # Member variables for full paths
    data_dir: str = field(init=False)
    config_path: str = field(init=False)
    sequences_fa: str = field(init=False)

    def __post_init__(self):
        self.data_dir = os.path.join(self.db_dir, "data", self.db_name)
        self.config_path = os.path.join(self.db_dir, self.config_filename)
        self.sequences_fa = os.path.join(self.data_dir, self.sequences_fa_filename)

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> "CustomDatabase":
        return cls(
            db_dir=args.db_dir,
            codon_table=getattr(args, "codon_table", DEFAULT_CODON_TABLE),
        )

    def find_annotation_file(self) -> Optional[str]:
        """Locate genes.gbk, genes.gff, or genes.gtf in data_dir."""
        for candidate in [self.genes_gbk_filename, self.genes_gff_filename, self.genes_gtf_filename]:
            path = os.path.join(self.data_dir, candidate)
            if os.path.exists(path):
                return path
        return None


@dataclass
class DbBuildConfig:
    """Parameters for building a custom database."""
    ref: str
    annotation_file: str
    db: CustomDatabase
    env: EnvironmentConfig
    chr_map: Optional[str] = None

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> "DbBuildConfig":
        return cls(
            ref=args.ref,
            annotation_file=args.annotation_file,
            db=CustomDatabase.from_args(args),
            env=EnvironmentConfig.from_args(args),
            chr_map=getattr(args, "chr_map", None),
        )


@dataclass
class AnnotationConfig:
    """Parameters for running batch SnpEff annotation."""
    input_dir: str
    strain_list: str
    output_dir: str
    db: CustomDatabase
    env: EnvironmentConfig
    ud: int = DEFAULT_UD

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> "AnnotationConfig":
        return cls(
            input_dir=args.input_dir,
            strain_list=args.strain_list,
            output_dir=args.output_dir,
            db=CustomDatabase.from_args(args),
            env=EnvironmentConfig.from_args(args),
            ud=getattr(args, "ud", DEFAULT_UD),
        )

    def get_strains(self) -> List[str]:
        """Read and return list of strains from strain_list file."""
        if not os.path.exists(self.strain_list):
            raise FileNotFoundError(f"Strain list file not found: {self.strain_list}")
        with open(self.strain_list, "r") as f:
            return [line.strip() for line in f if line.strip()]


# ==============================================================================
# Command Line Parser
# ==============================================================================
def parse_args(sys_args: Optional[List[str]] = None) -> argparse.Namespace:
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("--snpeff_env", default="snpeff", help="Conda environment for SnpEff")
    parent_parser.add_argument("--bcftools_env", default="snippy", help="Conda environment for bcftools")
    parent_parser.add_argument("--snippy_env", default="snippy", help="Conda environment for Snippy (e.g., snippy-vcf_to_tab)")
    parent_parser.add_argument("--bioperl_env", default="perl-bioperl-core", help="Conda environment for BioPerl tools (e.g., bp_genbank2gff3.pl)")

    parser = argparse.ArgumentParser(description="Utilities & Pipeline for SnpEff: Build database and run batch annotation.")
    subparsers = parser.add_subparsers(dest="command", required=True, help="Sub-commands")

    # 1. run subcommand (Full pipeline: DB Build + Strain Annotation)
    run_pipeline_parser = subparsers.add_parser(
        "run", parents=[parent_parser], 
        help="Full pipeline: Build a custom SnpEff database AND run batch annotation on strains"
    )
    run_pipeline_parser.add_argument("ref", help="Path to reference FASTA genome file")
    run_pipeline_parser.add_argument("annotation_file", help="Path to the custom GFF/GTF/Genbank annotation file")
    run_pipeline_parser.add_argument("input_dir", help="Directory containing strain subdirectories with snps.vcf")
    run_pipeline_parser.add_argument("strain_list", help="Text file containing one strain name per line")
    run_pipeline_parser.add_argument("output_dir", help="Directory to save SnpEff output")
    run_pipeline_parser.add_argument("--db_dir", default=DEFAULT_DB_DIR, help=f"Folder to store/read custom DB and config (default: {DEFAULT_DB_DIR})")
    run_pipeline_parser.add_argument("--codon_table", default=DEFAULT_CODON_TABLE, help=f"Codon table (e.g., {DEFAULT_CODON_TABLE})")
    run_pipeline_parser.add_argument("--chr_map", help="Tab-separated file for chromosome mapping (annotation_name ref_name)")
    run_pipeline_parser.add_argument("--ud", "-ud", type=int, default=DEFAULT_UD, help=f"Upstream/downstream interval length in bases (default: {DEFAULT_UD})")

    # 2. annotate subcommand (Batch annotation using existing DB)
    ann_parser = subparsers.add_parser("annotate", parents=[parent_parser], help="Run batch SnpEff annotation on strains using existing DB")
    ann_parser.add_argument("input_dir", help="Directory containing strain subdirectories with snps.vcf")
    ann_parser.add_argument("strain_list", help="Text file containing one strain name per line")
    ann_parser.add_argument("output_dir", help="Directory to save SnpEff output")
    ann_parser.add_argument("--db_dir", default=DEFAULT_DB_DIR, help=f"Folder to store/read custom DB and config (default: {DEFAULT_DB_DIR})")
    ann_parser.add_argument("--ud", "-ud", type=int, default=DEFAULT_UD, help=f"Upstream/downstream interval length in bases (default: {DEFAULT_UD})")

    # 3. db-build subcommand (Build custom SnpEff DB only)
    build_parser = subparsers.add_parser("db-build", parents=[parent_parser], help="Build a custom SnpEff database only")
    build_parser.add_argument("ref", help="Path to reference FASTA genome file")
    build_parser.add_argument("annotation_file", help="Path to the custom GFF/GTF/Genbank annotation file")
    build_parser.add_argument("--db_dir", default=DEFAULT_DB_DIR, help=f"Folder to store the custom DB and config (default: {DEFAULT_DB_DIR})")
    build_parser.add_argument("--codon_table", default=DEFAULT_CODON_TABLE, help=f"Codon table (e.g., {DEFAULT_CODON_TABLE})")
    build_parser.add_argument("--chr_map", help="Tab-separated file for chromosome mapping (annotation_name ref_name)")

    return parser.parse_args(sys_args)


# ==============================================================================
# Helper Functions
# ==============================================================================
def run_command(
    cmd: List[str], 
    env_name: str = "snippy", 
    check: bool = True, 
    stdout: Optional[IO] = None
) -> subprocess.CompletedProcess:
    """Run a command using conda run in the specified environment."""
    env_flag = "-p" if (os.path.isabs(env_name) or env_name.startswith("./") or env_name.startswith("../")) else "-n"
    full_cmd = ["conda", "run", "--no-capture-output", env_flag, env_name] + cmd

    if not stdout:
        logger.info(f"Running: {' '.join(full_cmd)}")
        result = subprocess.run(full_cmd, capture_output=False, text=True)
    else:
        logger.info(f"Running: {' '.join(full_cmd)} > {stdout.name if hasattr(stdout, 'name') else 'PIPE'}")
        result = subprocess.run(full_cmd, stdout=stdout)

    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed (exit code {result.returncode}): {' '.join(full_cmd)}")
    return result


def get_base_config_path(env_name: str) -> Optional[str]:
    """Locate the base snpEff.config template from the specified conda environment."""
    try:
        result = subprocess.run(
            ["conda", "run", "--no-capture-output", "-n", env_name, "which", "snpEff"],
            capture_output=True, text=True, check=True
        )
        snpeff_bin = result.stdout.strip()
        env_root = os.path.dirname(os.path.dirname(snpeff_bin))
        share_dir = os.path.join(env_root, "share")
        if os.path.exists(share_dir):
            for entry in os.listdir(share_dir):
                if entry.startswith("snpeff") and os.path.isdir(os.path.join(share_dir, entry)):
                    candidate = os.path.join(share_dir, entry, CONFIG_FILE_NAME)
                    if os.path.exists(candidate):
                        return candidate
    except Exception as e:
        logger.warning(f"Could not automatically locate base {CONFIG_FILE_NAME} in '{env_name}': {e}")
    return None


def get_annotation_chromosomes(annotation_file: str) -> List[str]:
    """Extract chromosome / VERSION / LOCUS names from GenBank, GFF, or GTF files."""
    if not annotation_file or not os.path.exists(annotation_file):
        return []

    chromosomes: List[str] = []
    ext = os.path.splitext(annotation_file)[1].lower()
    is_gbk = ext in [".gbk", ".gb", ".genbank", ".gbff"]

    try:
        with open(annotation_file, "r") as f:
            if is_gbk:
                current_locus = None
                current_version = None
                for line in f:
                    if line.startswith("LOCUS"):
                        if current_version and current_version not in chromosomes:
                            chromosomes.append(current_version)
                        elif current_locus and current_locus not in chromosomes:
                            chromosomes.append(current_locus)
                        parts = line.split()
                        current_locus = parts[1] if len(parts) > 1 else None
                        current_version = None
                    elif line.startswith("VERSION"):
                        parts = line.split()
                        if len(parts) > 1 and parts[1] != ".":
                            current_version = parts[1]
                if current_version and current_version not in chromosomes:
                    chromosomes.append(current_version)
                elif current_locus and current_locus not in chromosomes:
                    chromosomes.append(current_locus)
            else:
                for line in f:
                    if line.startswith("##sequence-region"):
                        parts = line.split()
                        if len(parts) > 1 and parts[1] not in chromosomes:
                            chromosomes.append(parts[1])
                    elif not line.startswith("#"):
                        parts = line.split("\t")
                        if len(parts) > 0 and parts[0] and parts[0] not in chromosomes:
                            chromosomes.append(parts[0])
    except Exception as e:
        logger.warning(f"Could not parse annotation file for chromosome names: {e}")

    return chromosomes


def get_fasta_chromosomes(fasta_path: str) -> List[str]:
    """Extract chromosome names (first word of > header) from FASTA file."""
    if not fasta_path or not os.path.exists(fasta_path):
        return []
    chromosomes: List[str] = []
    try:
        with open(fasta_path, "r") as f:
            for line in f:
                if line.startswith(">"):
                    cid = line[1:].strip().split()[0]
                    if cid and cid not in chromosomes:
                        chromosomes.append(cid)
    except Exception as e:
        logger.warning(f"Could not parse FASTA file for chromosome names: {e}")
    return chromosomes


def parse_chr_map(
    chr_map_path: str, 
    ann_chromosomes: List[str], 
    ref_chromosomes: List[str]
) -> Dict[str, str]:
    """
    Parse a tab-separated chromosome map file.
    Expected format per line: [annotation_name] [ref_fasta_name]
    Raises ValueError if annotation name does not exist in annotation file.
    """
    if not chr_map_path or not os.path.exists(chr_map_path):
        raise FileNotFoundError(f"Chromosome map file not found: {chr_map_path}")

    mapping: Dict[str, str] = {}
    with open(chr_map_path, "r") as f:
        for line_no, raw_line in enumerate(f, 1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 2:
                raise ValueError(
                    f"Invalid format in {chr_map_path} line {line_no}: expected '[annotation_name] [ref_name]', got '{line}'"
                )
            ann_name, ref_name = parts[0], parts[1]

            # Rule1: Validate that ann_name is present in annotation chromosomes (exact match)
            if ann_name not in ann_chromosomes:
                raise ValueError(
                    f"Chromosome '{ann_name}' specified in {chr_map_path} (line {line_no}) "
                    f"not found in annotation file chromosomes {ann_chromosomes}. "
                    f"Format must be: [annotation_name] [ref_name]"
                )

            # Rule 2: Duplicate source chromosome (key duplication)
            if ann_name in mapping:
                raise ValueError(
                    f"Duplicate mapping in {chr_map_path} (line {line_no}): annotation chromosome '{ann_name}' "
                    f"is specified more than once."
                )

            # Rule 3: Validate that ref_name is present in reference FASTA chromosomes
            if ref_name not in ref_chromosomes:
                raise ValueError(
                    f"Target chromosome '{ref_name}' specified in {chr_map_path} (line {line_no}) "
                    f"not found in reference FASTA chromosomes {ref_chromosomes}."
                )

            # Rule 4: Multiple source chromosomes mapped to the same target chromosome (value collision)
            if ref_name in mapping.values():
                raise ValueError(
                    f"Duplicate target mapping in {chr_map_path} (line {line_no}): target chromosome '{ref_name}' "
                    f"is already mapped from another chromosome."
                )

            mapping[ann_name] = ref_name

    if not mapping:
        raise ValueError(f"No valid chromosome mappings found in {chr_map_path}")

    return mapping


def adapt_annotation_for_db(annotation_file: str, output_file: str, mapping: Dict[str, str]) -> None:
    """
    Adapt chromosome names in annotation file (GBK, GFF, GTF) according to mapping (annotation_name -> ref_name).
    """
    if not mapping:
        shutil.copy(annotation_file, output_file)
        return

    ext = os.path.splitext(annotation_file)[1].lower()
    is_gbk = ext in [".gbk", ".gb", ".genbank", ".gbff"]

    if is_gbk:
        def _write_gbk_record(lines: List[str], f_out: IO[str], mapping: Dict[str, str]) -> None:
            locus_val = None
            version_val = None
            for l in lines:
                if l.startswith("LOCUS"):
                    parts = l.split()
                    if len(parts) > 1:
                        locus_val = parts[1]
                elif l.startswith("VERSION"):
                    parts = l.split()
                    if len(parts) > 1 and parts[1] != ".":
                        version_val = parts[1]

            target = None
            if version_val and version_val in mapping:
                target = mapping[version_val]
            elif locus_val and locus_val in mapping:
                target = mapping[locus_val]


            for l in lines:
                if target:
                    if l.startswith("LOCUS"):
                        parts = l.split()
                        prefix = "LOCUS       "
                        trailing = " ".join(parts[2:]) if len(parts) > 2 else ""
                        l = f"{prefix}{target:<20} {trailing}\n"
                    elif l.startswith("VERSION"):
                        l = f"VERSION     {target}\n"
                f_out.write(l)

        with open(annotation_file, "r") as f_in, open(output_file, "w") as f_out:
            rec_lines = []
            for line in f_in:
                rec_lines.append(line)
                if line.startswith("//"):
                    _write_gbk_record(rec_lines, f_out, mapping)
                    rec_lines = []
            if rec_lines:
                _write_gbk_record(rec_lines, f_out, mapping)

    else:
        # GFF3 or GTF
        with open(annotation_file, "r") as f_in, open(output_file, "w") as f_out:
            for line in f_in:
                if line.startswith("##sequence-region"):
                    parts = line.strip().split()
                    if len(parts) > 1 and parts[1] in mapping:
                        parts[1] = mapping[parts[1]]
                        line = " ".join(parts) + "\n"
                elif not line.startswith("#"):
                    parts = line.split("\t")
                    if len(parts) >= 9:
                        seq_id = parts[0]
                        if seq_id in mapping:
                            parts[0] = mapping[seq_id]
                            line = "\t".join(parts)
                f_out.write(line)


def get_vcf_chromosomes(vcf_file: str) -> List[str]:
    """Extract chromosome names from VCF header or #CHROM records."""
    if not vcf_file or not os.path.exists(vcf_file):
        return []

    chromosomes: List[str] = []
    try:
        with open(vcf_file, "r") as f:
            for line in f:
                if line.startswith("##contig="):
                    match = re.search(r"ID=([^,>]+)", line)
                    if match:
                        cid = match.group(1)
                        if cid not in chromosomes:
                            chromosomes.append(cid)
                elif not line.startswith("#"):
                    parts = line.split("\t")
                    if len(parts) > 0 and parts[0] and parts[0] not in chromosomes:
                        chromosomes.append(parts[0])
    except Exception as e:
        logger.warning(f"Could not parse VCF file for chromosome names: {e}")

    return chromosomes


def rename_chromosomes_bcftools(
    input_vcf: str, 
    output_vcf: str, 
    chr_map_path: str, 
    env_name: str = "snippy"
) -> None:
    """Rename chromosomes in VCF using bcftools annotate --rename-chrs."""
    if not os.path.exists(chr_map_path):
        logger.warning(f"Chromosome map file not found: {chr_map_path}")
        shutil.copy(input_vcf, output_vcf)
        return

    cmd = ["bcftools", "annotate", "--no-version", "--rename-chrs", chr_map_path, input_vcf, "-o", output_vcf]
    run_command(cmd, env_name=env_name, check=True)


def generate_snippy_tab(
    vcf_path: str,
    output_tab_path: str,
    annotation_file: str,
    ref_fa: str,
    env_name: str = "snippy",
    bioperl_env: str = "perl-bioperl-core"
) -> None:
    """Generate Snippy-style snps.tab table using bp_genbank2gff3 and snippy-vcf_to_tab."""
    if not vcf_path or not os.path.exists(vcf_path):
        return

    ext = os.path.splitext(annotation_file)[1].lower() if annotation_file else ""
    is_gbk = ext in [".gbk", ".gb", ".genbank", ".gbff"]

    # Use a persistent debug directory in output folder for verification
    out_dir = os.path.dirname(output_tab_path)
    debug_dir = os.path.join(out_dir, "debug_gff_files")
    os.makedirs(debug_dir, exist_ok=True)

    ref_gff = annotation_file if not is_gbk else None
    target_vcf = vcf_path

    if is_gbk and annotation_file and os.path.exists(annotation_file):
        # Convert GBK to GFF3 using bp_genbank2gff3.pl (perl-bioperl-core) or bp_genbank2gff3
        #    -noinfer: don't infer exon/mRNA subfeatures (L1022: passed to unflatten_seq)
        #    --filter gene: remove gene features from output (L1049: sub filter)
        #    --filter exon: remove exon features from output
        derived_gff = os.path.join(debug_dir, "derived_ref.gff")
        bp_cmd = ["perl", "-S", "bp_genbank2gff3.pl", "-noinfer", "--filter", "gene", "--filter", "exon", "-out", "stdout", annotation_file]
        try:
            with open(derived_gff, "w") as out_gff_f:
                result = run_command(bp_cmd, env_name=bioperl_env, stdout=out_gff_f, check=False)
                if result.returncode != 0:
                    fallback_cmd = ["perl", "-S", "bp_genbank2gff3", "-noinfer", "--filter", "gene", "--filter", "exon", "-out", "stdout", annotation_file]
                    out_gff_f.seek(0)
                    out_gff_f.truncate(0)
                    run_command(fallback_cmd, env_name=bioperl_env, stdout=out_gff_f, check=False)
            if os.path.exists(derived_gff) and os.path.getsize(derived_gff) > 0:
                ref_gff = derived_gff
        except Exception as e:
            logger.warning(f"Could not convert GBK to GFF3: {e}")

    # Clean FASTA/GFF/VCF:
    # 1. Filter out GFF features that do not exist in the reference FASTA file
    # 2. Auto-map VCF CHROM names to reference FASTA sequence IDs if they differ.
    valid_seq_ids = []
    if ref_fa and os.path.exists(ref_fa):
        clean_fa = os.path.join(debug_dir, "clean_ref.fa")
        try:
            with open(ref_fa, "r") as f_in, open(clean_fa, "w") as f_out:
                for line in f_in:
                    if line.startswith(">"):
                        seq_id = line.split()[0][1:]
                        if seq_id not in valid_seq_ids:
                            valid_seq_ids.append(seq_id)
                        f_out.write(f">{seq_id}\n")
                    else:
                        f_out.write(line)
            ref_fa = clean_fa
        except Exception as e:
            logger.warning(f"Could not prepare cleaned FASTA for Snippy tab: {e}")

    valid_seq_set = set(valid_seq_ids)

    if ref_gff and os.path.exists(ref_gff):
        clean_gff = os.path.join(debug_dir, "clean_ref.gff")
        try:
            with open(ref_gff, "r") as f_in, open(clean_gff, "w") as f_out:
                for line in f_in:
                    if line.startswith("#"):
                        f_out.write(line)
                        continue
                    parts = line.split("\t")
                    if len(parts) < 9:
                        f_out.write(line)
                        continue
                    seq_id = parts[0]
                    ftype = parts[2]
                    # Filter out sequence IDs absent in FASTA
                    if valid_seq_set and seq_id not in valid_seq_set:
                        continue
                    
                    # Parse attributes supporting both GFF3 (key=value) and GTF (key "value")
                    raw_attrs = parts[8].strip().split(";")
                    attr_dict = {}
                    attr_order = []
                    for item in raw_attrs:
                        item = item.strip()
                        if not item:
                            continue
                        if "=" in item:
                            k, v = item.split("=", 1)
                            k = k.strip()
                            v = v.strip()
                            attr_dict[k] = v
                            if k not in attr_order:
                                attr_order.append(k)
                        else:
                            m = re.match(r'^(\S+)\s+(.+)$', item)
                            if m:
                                k = m.group(1).strip()
                                val_part = m.group(2).strip()
                                vals = re.findall(r'"([^"]*)"', val_part)
                                v = ",".join(vals) if vals else val_part.strip('"')
                                attr_dict[k] = v
                                if k not in attr_order:
                                    attr_order.append(k)

                    # For CDS features, ensure 'gene' and 'locus_tag' attributes exist
                    if ftype.upper() == "CDS":
                        if "gene" not in attr_dict or not attr_dict["gene"]:
                            for candidate in ["Name", "gene_name", "gene_id"]:
                                if candidate in attr_dict and attr_dict[candidate]:
                                    attr_dict["gene"] = attr_dict[candidate]
                                    if "gene" not in attr_order:
                                        attr_order.append("gene")
                                    break
                        if "locus_tag" not in attr_dict or not attr_dict["locus_tag"]:
                            for candidate in ["gene_id", "transcript_id", "Name"]:
                                if candidate in attr_dict and attr_dict[candidate]:
                                    attr_dict["locus_tag"] = attr_dict[candidate]
                                    if "locus_tag" not in attr_order:
                                        attr_order.append("locus_tag")
                                    break

                    # Reconstruct attributes string in GFF3 format (key=value) so snippy-vcf_to_tab can parse it
                    if attr_order:
                        new_attr_str = ";".join([f"{k}={attr_dict[k]}" for k in attr_order])
                        parts[8] = new_attr_str + "\n"
                        line = "\t".join(parts)

                    f_out.write(line)
            ref_gff = clean_gff
        except Exception as e:
            logger.warning(f"Could not prepare cleaned GFF for Snippy tab: {e}")

    if vcf_path and os.path.exists(vcf_path):
        clean_vcf = os.path.join(debug_dir, "clean_var.vcf")
        
        # Determine target reference CHROM name for single-chromosome genomes
        target_ref_chr = valid_seq_ids[0] if len(valid_seq_ids) == 1 else None

        try:
            with open(vcf_path, "r") as f_in, open(clean_vcf, "w") as f_out:
                for line in f_in:
                    if line.startswith("#"):
                        f_out.write(line)
                        continue
                    parts = line.split("\t")
                    if len(parts) < 8:
                        f_out.write(line)
                        continue
                    vcf_chr = parts[0]
                    # Remap CHROM to reference ID if single chromosome and name differs
                    if target_ref_chr and vcf_chr not in valid_seq_set:
                        parts[0] = target_ref_chr
                    f_out.write("\t".join(parts))
            target_vcf = clean_vcf
        except Exception as e:
            logger.warning(f"Could not prepare cleaned VCF for Snippy tab: {e}")

    # Check if required FASTA and GFF3 reference files were successfully prepared
    if not ref_fa or not ref_gff:
        strain_name = os.path.basename(output_tab_path).replace(".snps.tab", "")
        logger.warning(
            f"Could not prepare FASTA/GFF3 for Snippy tab table for strain '{strain_name}'. Skipping."
        )
        return

    # Use perl -X to suppress uninitialized value warnings from snippy-vcf_to_tab
    cmd = ["perl", "-X", "-S", "snippy-vcf_to_tab"]
    if ref_fa:
        cmd += ["--ref", ref_fa]
    if ref_gff:
        cmd += ["--gff", ref_gff]
    cmd += ["--vcf", target_vcf]

    try:
        with open(output_tab_path, "w") as out_f:
            run_command(cmd, env_name=env_name, check=False, stdout=out_f)
        if os.path.exists(output_tab_path) and os.path.getsize(output_tab_path) > 0:
            logger.info(f"Successfully created Snippy tab table: {output_tab_path}")
    except Exception as e:
        logger.warning(f"Could not generate Snippy tab table: {e}")


# ==============================================================================
# Pipeline Execution Functions
# ==============================================================================
def build_custom_db(config: DbBuildConfig) -> None:
    """Builds a custom SnpEff database from reference FASTA and annotation files."""
    db = config.db
    env = config.env

    logger.info(f"=== Building Custom Database: {db.db_name} ===")
    os.makedirs(db.data_dir, exist_ok=True)
    
    # 1. Check Fasta and Annotation
    if not os.path.exists(config.ref):
        raise FileNotFoundError(f"Reference FASTA file not found: {config.ref}")

    if not os.path.exists(config.annotation_file):
        raise FileNotFoundError(f"Annotation file not found: {config.annotation_file}")

    ref_chromosomes = get_fasta_chromosomes(config.ref)
    ann_chromosomes = get_annotation_chromosomes(config.annotation_file)
    logger.info(f"Reference FASTA chromosomes: {ref_chromosomes}")
    logger.info(f"Annotation chromosomes: {ann_chromosomes}")

    if not ref_chromosomes:
        raise ValueError(f"No chromosomes found in reference FASTA file: {config.ref}")
    if not ann_chromosomes:
        raise ValueError(f"No chromosomes found in annotation file: {config.annotation_file}")

    # Determine mapping (annotation_name -> ref_fasta_name)
    mapping: Dict[str, str] = {}
    if config.chr_map:
        mapping = parse_chr_map(config.chr_map, ann_chromosomes, ref_chromosomes)
        logger.info(f"Using explicit chromosome mapping for DB build: {mapping}")
    elif len(ref_chromosomes) == 1 and len(ann_chromosomes) == 1:
        if ref_chromosomes[0] != ann_chromosomes[0]:
            mapping = {ann_chromosomes[0]: ref_chromosomes[0]}
            logger.info(f"Auto-mapped single chromosome for DB build: '{ann_chromosomes[0]}' -> '{ref_chromosomes[0]}'")
    else:
        common = set(ref_chromosomes).intersection(set(ann_chromosomes))
        if not common:
            raise ValueError(
                f"No matching chromosome names between Reference FASTA {ref_chromosomes} and Annotation {ann_chromosomes}. "
                f"Please specify --chr_map with '[annotation_name] [ref_name]'."
            )

    ext = os.path.splitext(config.annotation_file)[1].lower()
    is_gtf = ext == ".gtf"
    is_gbk = ext in [".gbk", ".gb", ".genbank", ".gbff"]

    if is_gbk:
        ann_target = os.path.join(db.data_dir, db.genes_gbk_filename)
    elif is_gtf:
        ann_target = os.path.join(db.data_dir, db.genes_gtf_filename)
    else:
        ann_target = os.path.join(db.data_dir, db.genes_gff_filename)
        
    ref_target = db.sequences_fa

    adapt_annotation_for_db(config.annotation_file, ann_target, mapping)
    shutil.copy(config.ref, ref_target)
    logger.info(f"Prepared database files in {db.data_dir}")

    # 2. Setup Config
    local_config = db.config_path
    base_config = get_base_config_path(env.snpeff_env)
    if base_config and os.path.exists(base_config):
        shutil.copy(base_config, local_config)
        logger.info(f"Copied base config from {base_config}")
    else:
        with open(local_config, "w") as f:
            f.write("# SnpEff configuration\n")
        logger.warning(f"Base config not found. Created minimal config at {local_config}")

    with open(local_config, "a") as f:
        f.write(f"\n# Custom Genome Entry\ndata.dir = ./data\n{db.db_name}.genome : {db.db_name} {db.codon_table}\n")
    logger.info(f"Configured database in {local_config}")

    # 3. Run Build
    if is_gbk:
        build_flag = "-genbank"
    elif is_gtf:
        build_flag = "-gtf22"
    else:
        build_flag = "-gff3"

    build_cmd = [
        "snpEff", "build", 
        "-c", local_config, 
        build_flag, 
        "-noCheckCds", "-noCheckProtein", 
        "-v", db.db_name
    ]
    
    run_command(build_cmd, env_name=env.snpeff_env)
    logger.info(f"Database build processed for {db.db_name}")


def run_annotation(config: AnnotationConfig) -> None:
    """Executes batch SnpEff annotation and Snippy tab table generation for strains."""
    db = config.db
    env = config.env

    if not os.path.exists(db.db_dir):
        raise FileNotFoundError(f"Database directory not found: {db.db_dir}. Please run db-build first.")

    if not os.path.exists(db.data_dir):
        raise FileNotFoundError(f"Database directory not found: {db.data_dir}. Please run db-build first.")

    ref_fa = db.sequences_fa
    if not os.path.exists(ref_fa):
        raise FileNotFoundError(f"Reference FASTA not found in database: {ref_fa}")

    ann_file = db.find_annotation_file()
    if not ann_file:
        raise FileNotFoundError(
            f"Annotation file ({db.genes_gbk_filename}/{db.genes_gff_filename}/{db.genes_gtf_filename}) "
            f"not found in database: {db.data_dir}"
        )

    config_path = db.config_path
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found in database directory: {config_path}. Please run db-build first.")

    if not os.path.exists(config.output_dir):
        os.makedirs(config.output_dir, exist_ok=True)
        logger.info(f"Created output directory: {config.output_dir}")

    strains = config.get_strains()
    ann_chromosomes = get_annotation_chromosomes(ann_file)

    logger.info(f"Processing {len(strains)} strains using SnpEff DB: {db.db_name}")
    logger.info(f"Using config: {config_path}")
    if ann_chromosomes:
        logger.info(f"Reference annotation chromosomes: {ann_chromosomes}")

    for strain in strains:
        _process_single_strain(strain, config, ann_file, ref_fa, config_path, ann_chromosomes)

    logger.info("Batch processing completed.")


def _process_single_strain(
    strain: str, 
    config: AnnotationConfig, 
    ann_file: str, 
    ref_fa: str, 
    config_path: str, 
    ann_chromosomes: List[str]
) -> None:
    """Annotate and generate Snippy tab table for a single strain."""
    db = config.db
    env = config.env

    logger.info(f"--- Processing strain: {strain} ---")
    input_vcf = os.path.join(config.input_dir, strain, "snps.vcf")
    output_vcf = os.path.join(config.output_dir, f"{strain}.ann.vcf")
    stats_file = os.path.join(config.output_dir, f"{strain}_summary.csv")
    output_tab = os.path.join(config.output_dir, f"{strain}.snps.tab")

    if not os.path.exists(input_vcf):
        logger.warning(f"Input VCF not found for {strain}: {input_vcf}")
        return

    vcf_chromosomes = get_vcf_chromosomes(input_vcf)
    if vcf_chromosomes and ann_chromosomes:
        common_chromosomes = set(ann_chromosomes).intersection(set(vcf_chromosomes))
        if not common_chromosomes:
            logger.warning(
                f"No matching chromosome names between VCF {vcf_chromosomes} and Database {ann_chromosomes} "
                f"for strain '{strain}'. Skipping strain annotation."
            )
            return

    # SnpEff command (VCF is used directly since DB was built matching reference FASTA)
    snpeff_cmd = [
        "snpEff",
        "-c", config_path,
        "-v", 
        "-ud", str(config.ud),
        "-i", "vcf", 
        "-o", "vcf", 
        "-csvStats", stats_file,
        db.db_name, 
        input_vcf
    ]

    with open(output_vcf, "w") as out_f:
        run_command(snpeff_cmd, env_name=env.snpeff_env, stdout=out_f)
    logger.info(f"Successfully created: {output_vcf}")
    logger.info(f"Summary report: {stats_file}")

    # Clean up default html report generated in current working directory
    if os.path.exists("snpEff_summary.html"):
        try:
            os.remove("snpEff_summary.html")
        except OSError:
            pass

    # Generate Snippy snps.tab table using snippy_env and bioperl_env
    generate_snippy_tab(
        output_vcf, output_tab, ann_file, ref_fa,
        env_name=env.snippy_env,
        bioperl_env=env.bioperl_env
    )


# ==============================================================================
# Subcommand Handlers & Main
# ==============================================================================
def handle_run(build_config: DbBuildConfig, ann_config: AnnotationConfig) -> None:
    logger.info("=========================================")
    logger.info(" Starting Full Pipeline (Build + Annotate)")
    logger.info("=========================================")

    build_custom_db(build_config)
    run_annotation(ann_config)


def main() -> None:
    try:
        args = parse_args()
        if args.command == "db-build":
            build_config = DbBuildConfig.from_args(args)
            build_custom_db(build_config)
        elif args.command == "annotate":
            ann_config = AnnotationConfig.from_args(args)
            run_annotation(ann_config)
        elif args.command == "run":
            build_config = DbBuildConfig.from_args(args)
            ann_config = AnnotationConfig.from_args(args)
            handle_run(build_config, ann_config)
    except Exception as e:
        logger.error(f"{e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
