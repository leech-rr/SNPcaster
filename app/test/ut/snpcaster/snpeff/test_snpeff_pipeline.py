#!/usr/bin/env python3
import unittest
import os
import sys
import subprocess
import shutil
import tempfile
import datetime
from typing import Optional, Dict, List
from unittest.mock import patch

import snpeff_pipeline
from snpeff_pipeline import (
    CustomDatabase,
    DbBuildConfig,
    AnnotationConfig,
    CUSTOM_DB_NAME,
    CONFIG_FILE_NAME,
    DEFAULT_UD,
)


class TestParseArgs(unittest.TestCase):
    def test_run_subcommand(self):
        cmd = [
            "run", "ref.fasta", "ref.gbk", "./input_dir", "strains.txt", "./output_dir",
            "--db_dir", "./my_db_dir", "--snpeff_env", "my_snpeff_env", "--codon_table", "[Bacteria]"
        ]
        args = snpeff_pipeline.parse_args(cmd)
        self.assertEqual(args.command, "run")
        self.assertEqual(args.ref, "ref.fasta")
        self.assertEqual(args.annotation_file, "ref.gbk")
        self.assertEqual(args.input_dir, "./input_dir")
        self.assertEqual(args.strain_list, "strains.txt")
        self.assertEqual(args.output_dir, "./output_dir")
        self.assertEqual(args.db_dir, "./my_db_dir")
        self.assertEqual(args.snpeff_env, "my_snpeff_env")
        self.assertEqual(args.codon_table, "[Bacteria]")
        self.assertEqual(args.ud, DEFAULT_UD)

    def test_annotate_subcommand(self):
        cmd = ["annotate", "./input_dir", "strains.txt", "./output_dir", "--db_dir", "./my_db_dir"]
        args = snpeff_pipeline.parse_args(cmd)
        self.assertEqual(args.command, "annotate")
        self.assertEqual(args.input_dir, "./input_dir")
        self.assertEqual(args.strain_list, "strains.txt")
        self.assertEqual(args.output_dir, "./output_dir")
        self.assertEqual(args.db_dir, "./my_db_dir")
        self.assertEqual(args.snpeff_env, "snpeff")
        self.assertEqual(args.ud, DEFAULT_UD)

    def test_ud_argument(self):
        args1 = snpeff_pipeline.parse_args(["run", "ref.fa", "ref.gbk", "in", "strains.txt", "out", "--ud", "500"])
        self.assertEqual(args1.ud, 500)
        args2 = snpeff_pipeline.parse_args(["annotate", "in", "strains.txt", "out", "-ud", "100"])
        self.assertEqual(args2.ud, 100)
        config1 = AnnotationConfig.from_args(args1)
        self.assertEqual(config1.ud, 500)
        config2 = AnnotationConfig.from_args(args2)
        self.assertEqual(config2.ud, 100)

    def test_db_build_subcommand(self):
        cmd = ["db-build", "ref.fasta", "ref.gbk", "--db_dir", "./out_db", "--codon_table", "[Bacteria]"]
        args = snpeff_pipeline.parse_args(cmd)
        self.assertEqual(args.command, "db-build")
        self.assertEqual(args.ref, "ref.fasta")
        self.assertEqual(args.annotation_file, "ref.gbk")
        self.assertEqual(args.db_dir, "./out_db")
        self.assertEqual(args.codon_table, "[Bacteria]")


class TestDataClasses(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_custom_database_paths(self):
        db = CustomDatabase(db_dir=self.temp_dir, codon_table="[Bacteria]")
        expected_data_dir = os.path.join(self.temp_dir, "data", CUSTOM_DB_NAME)
        expected_config = os.path.join(self.temp_dir, CONFIG_FILE_NAME)
        expected_seq = os.path.join(expected_data_dir, "sequences.fa")

        self.assertEqual(db.data_dir, expected_data_dir)
        self.assertEqual(db.config_path, expected_config)
        self.assertEqual(db.sequences_fa, expected_seq)

    def test_find_annotation_file(self):
        db = CustomDatabase(db_dir=self.temp_dir, codon_table="[Bacteria]")
        os.makedirs(db.data_dir, exist_ok=True)
        gbk_path = os.path.join(db.data_dir, "genes.gbk")
        with open(gbk_path, "w") as f:
            f.write("LOCUS test 100 bp")

        self.assertEqual(db.find_annotation_file(), gbk_path)

    def test_annotation_config_get_strains(self):
        strain_file = os.path.join(self.temp_dir, "strains.txt")
        with open(strain_file, "w") as f:
            f.write("strain1\n\nstrain2\n")

        args = snpeff_pipeline.parse_args(["annotate", "./input", strain_file, "./output", "--db_dir", self.temp_dir])
        config = AnnotationConfig.from_args(args)
        strains = config.get_strains()
        self.assertEqual(strains, ["strain1", "strain2"])


class TestChromosomeParsing(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.gbk_file = os.path.join(self.temp_dir, "test.gbk")
        with open(self.gbk_file, "w") as f:
            f.write("LOCUS       contig_1             5000 bp    DNA\nVERSION     BA000007.3\n")

        self.fasta_file = os.path.join(self.temp_dir, "test.fasta")
        with open(self.fasta_file, "w") as f:
            f.write(">contig_fasta extra description\nATGC\n")

        self.vcf_file = os.path.join(self.temp_dir, "snps.vcf")
        with open(self.vcf_file, "w") as f:
            f.write("##fileformat=VCFv4.2\n##contig=<ID=O26_11368,length=5000>\n#CHROM\tPOS\tID\tREF\tALT\nO26_11368\t10\t.\tA\tT\n")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_get_annotation_chromosomes(self):
        chrs = snpeff_pipeline.get_annotation_chromosomes(self.gbk_file)
        self.assertEqual(chrs, ["BA000007.3"])

    def test_get_fasta_chromosomes(self):
        chrs = snpeff_pipeline.get_fasta_chromosomes(self.fasta_file)
        self.assertEqual(chrs, ["contig_fasta"])

    def test_get_vcf_chromosomes(self):
        chrs = snpeff_pipeline.get_vcf_chromosomes(self.vcf_file)
        self.assertEqual(chrs, ["O26_11368"])

    def test_parse_chr_map_valid(self):
        map_file = os.path.join(self.temp_dir, "map.txt")
        with open(map_file, "w") as f:
            f.write("BA000007.3\tchr_ref\n")
        mapping = snpeff_pipeline.parse_chr_map(map_file, ["BA000007.3"], ["chr_ref"])
        self.assertEqual(mapping.get("BA000007.3"), "chr_ref")

    def test_parse_chr_map_reverse_raises_error(self):
        map_file = os.path.join(self.temp_dir, "map_reverse.txt")
        with open(map_file, "w") as f:
            f.write("chr_ref\tBA000007.3\n")
        with self.assertRaises(ValueError):
            snpeff_pipeline.parse_chr_map(map_file, ["BA000007.3"], ["chr_ref"])

    def test_parse_chr_map_invalid_column_count(self):
        # 1 column
        map_file_1 = os.path.join(self.temp_dir, "map_1col.txt")
        with open(map_file_1, "w") as f:
            f.write("BA000007.3\n")
        with self.assertRaises(ValueError):
            snpeff_pipeline.parse_chr_map(map_file_1, ["BA000007.3"], ["chr_ref"])

        # 3 columns
        map_file_3 = os.path.join(self.temp_dir, "map_3col.txt")
        with open(map_file_3, "w") as f:
            f.write("BA000007.3\tchr_ref\textra_column\n")
        with self.assertRaises(ValueError):
            snpeff_pipeline.parse_chr_map(map_file_3, ["BA000007.3"], ["chr_ref"])

    def test_parse_chr_map_identical_names_allowed(self):
        """Source and target chromosome names are identical (allowed, explicit identity mapping)."""
        map_file = os.path.join(self.temp_dir, "map_identical.txt")
        with open(map_file, "w") as f:
            f.write("BA000007.3\tBA000007.3\n")
        mapping = snpeff_pipeline.parse_chr_map(map_file, ["BA000007.3"], ["BA000007.3"])
        self.assertEqual(mapping.get("BA000007.3"), "BA000007.3")

    def test_parse_chr_map_duplicate_source_raises_error(self):
        """Rule 2: Duplicate source (annotation) chromosome specified."""
        map_file = os.path.join(self.temp_dir, "map_dup_source.txt")
        with open(map_file, "w") as f:
            f.write("BA000007.3\tchr1\nBA000007.3\tchr2\n")
        with self.assertRaises(ValueError) as ctx:
            snpeff_pipeline.parse_chr_map(map_file, ["BA000007.3"], ["chr1", "chr2"])
        self.assertIn("is specified more than once", str(ctx.exception))

    def test_parse_chr_map_duplicate_target_raises_error(self):
        """Rule 3: Multiple annotation chromosomes mapped to the same reference chromosome."""
        map_file = os.path.join(self.temp_dir, "map_dup_target.txt")
        with open(map_file, "w") as f:
            f.write("BA000007.3\tchr1\nBA000008.1\tchr1\n")
        with self.assertRaises(ValueError) as ctx:
            snpeff_pipeline.parse_chr_map(map_file, ["BA000007.3", "BA000008.1"], ["chr1"])
        self.assertIn("is already mapped from another chromosome", str(ctx.exception))

    def test_parse_chr_map_target_not_in_ref_raises_error(self):
        """Rule 4: Target chromosome name does not exist in reference FASTA chromosomes."""
        map_file = os.path.join(self.temp_dir, "map_target_not_in_ref.txt")
        with open(map_file, "w") as f:
            f.write("BA000007.3\tunknown_ref\n")
        with self.assertRaises(ValueError) as ctx:
            snpeff_pipeline.parse_chr_map(map_file, ["BA000007.3"], ["chr_ref"])
        self.assertIn("not found in reference FASTA chromosomes", str(ctx.exception))

    def test_adapt_annotation_gbk(self):
        out_gbk = os.path.join(self.temp_dir, "adapted.gbk")
        snpeff_pipeline.adapt_annotation_for_db(self.gbk_file, out_gbk, {"BA000007.3": "chr_ref"})
        with open(out_gbk) as f:
            content = f.read()
        self.assertIn("VERSION     chr_ref", content)
        self.assertIn("LOCUS       chr_ref", content)

    def test_adapt_annotation_gff(self):
        gff_file = os.path.join(self.temp_dir, "test.gff")
        out_gff = os.path.join(self.temp_dir, "adapted.gff")
        with open(gff_file, "w") as f:
            f.write("##sequence-region BA000007.3 1 1000\nBA000007.3\tRefSeq\tgene\t1\t100\t.\t+\t.\tID=g1\n")
        snpeff_pipeline.adapt_annotation_for_db(gff_file, out_gff, {"BA000007.3": "chr_ref"})
        with open(out_gff) as f:
            content = f.read()
        self.assertIn("##sequence-region chr_ref 1 1000", content)
        self.assertTrue(content.splitlines()[1].startswith("chr_ref\tRefSeq"))


class TestBuildCustomDb(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.ref_fa = os.path.join(self.temp_dir, "test.fasta")
        with open(self.ref_fa, "w") as f:
            f.write(">contig_1\nATGC\n")

        self.annotation_file = os.path.join(self.temp_dir, "test.gbk")
        with open(self.annotation_file, "w") as f:
            f.write("LOCUS       test_locus 100 bp\nVERSION     test_locus.1\n")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    @patch("snpeff_pipeline.get_base_config_path", return_value=None)
    @patch("snpeff_pipeline.run_command")
    def test_build_custom_db(self, mock_run_command, mock_get_base_config):
        out_db_dir = os.path.join(self.temp_dir, "custom_db")
        args = snpeff_pipeline.parse_args([
            "db-build", self.ref_fa, self.annotation_file, 
            "--db_dir", out_db_dir
        ])
        config = DbBuildConfig.from_args(args)
        snpeff_pipeline.build_custom_db(config)
        mock_run_command.assert_called_once()
        self.assertTrue(os.path.exists(config.db.sequences_fa))
        self.assertTrue(os.path.exists(config.db.config_path))
        # Verify annotation was auto-adapted to contig_1
        with open(os.path.join(config.db.data_dir, "genes.gbk")) as f:
            content = f.read()
        self.assertIn("VERSION     contig_1", content)

    def test_build_custom_db_empty_fasta_raises_error(self):
        """Build custom DB with FASTA containing no chromosome headers raises ValueError."""
        empty_fa = os.path.join(self.temp_dir, "empty.fasta")
        with open(empty_fa, "w") as f:
            f.write("# no fasta headers here\n")
        out_db_dir = os.path.join(self.temp_dir, "custom_db_empty_fa")
        args = snpeff_pipeline.parse_args([
            "db-build", empty_fa, self.annotation_file,
            "--db_dir", out_db_dir
        ])
        config = DbBuildConfig.from_args(args)
        with self.assertRaises(ValueError) as ctx:
            snpeff_pipeline.build_custom_db(config)
        self.assertIn("No chromosomes found in reference FASTA file", str(ctx.exception))

    def test_build_custom_db_empty_annotation_raises_error(self):
        """Build custom DB with annotation containing no chromosome records raises ValueError."""
        empty_ann = os.path.join(self.temp_dir, "empty.gbk")
        with open(empty_ann, "w") as f:
            f.write("# empty gbk\n")
        out_db_dir = os.path.join(self.temp_dir, "custom_db_empty_ann")
        args = snpeff_pipeline.parse_args([
            "db-build", self.ref_fa, empty_ann,
            "--db_dir", out_db_dir
        ])
        config = DbBuildConfig.from_args(args)
        with self.assertRaises(ValueError) as ctx:
            snpeff_pipeline.build_custom_db(config)
        self.assertIn("No chromosomes found in annotation file", str(ctx.exception))


class TestRunAnnotation(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.input_dir = os.path.join(self.temp_dir, "input")
        self.output_dir = os.path.join(self.temp_dir, "output")
        self.strain_dir = os.path.join(self.input_dir, "strainA")
        os.makedirs(self.strain_dir)

        # Mock DB
        self.db_dir = os.path.join(self.temp_dir, "custom_db")
        self.data_dir = os.path.join(self.db_dir, "data", CUSTOM_DB_NAME)
        os.makedirs(self.data_dir)

        self.gbk_file = os.path.join(self.data_dir, "genes.gbk")
        with open(self.gbk_file, "w") as f:
            f.write("LOCUS       contig_1 100 bp\n")

        self.ref_fa = os.path.join(self.data_dir, "sequences.fa")
        with open(self.ref_fa, "w") as f:
            f.write(">contig_1\nATGC\n")

        self.config_path = os.path.join(self.db_dir, CONFIG_FILE_NAME)
        with open(self.config_path, "w") as f:
            f.write("# SnpEff configuration\n")

        # Mock VCF file
        self.vcf_file = os.path.join(self.strain_dir, "snps.vcf")
        with open(self.vcf_file, "w") as f:
            f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\ncontig_1\t10\t.\tA\tT\n")

        # Strain list
        self.strain_list = os.path.join(self.temp_dir, "strains.txt")
        with open(self.strain_list, "w") as f:
            f.write("strainA\n")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    @patch("snpeff_pipeline.generate_snippy_tab")
    @patch("snpeff_pipeline.run_command")
    def test_run_annotation_matching_contig(self, mock_run_command, mock_snippy_tab):
        """Test running annotation when VCF matches DB."""
        args = snpeff_pipeline.parse_args([
            "annotate", self.input_dir, self.strain_list, self.output_dir, 
            "--db_dir", self.db_dir
        ])
        config = AnnotationConfig.from_args(args)
        snpeff_pipeline.run_annotation(config)
        mock_run_command.assert_called_once()
        cmd = mock_run_command.call_args[0][0]
        self.assertIn("-ud", cmd)
        ud_idx = cmd.index("-ud")
        self.assertEqual(cmd[ud_idx + 1], "300")
        mock_snippy_tab.assert_called_once()

    @patch("snpeff_pipeline.run_command")
    def test_mismatched_contig_skip(self, mock_run_command):
        """Test skipping strain when VCF contig does not match DB."""
        # Update VCF to have non-matching contig
        with open(self.vcf_file, "w") as f:
            f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\nother_contig\t10\t.\tA\tT\n")

        args = snpeff_pipeline.parse_args([
            "annotate", self.input_dir, self.strain_list, self.output_dir, 
            "--db_dir", self.db_dir
        ])
        config = AnnotationConfig.from_args(args)
        snpeff_pipeline.run_annotation(config)
        mock_run_command.assert_not_called()


class TestGenerateSnippyTab(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.vcf_file = os.path.join(self.temp_dir, "test.vcf")
        with open(self.vcf_file, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
        self.out_tab = os.path.join(self.temp_dir, "test.snps.tab")
        self.ann_file = os.path.join(self.temp_dir, "genes.gbk")
        with open(self.ann_file, "w") as f:
            f.write("LOCUS test 100 bp\n")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_generate_snippy_tab_header_and_file_creation(self):
        """Verify tab file creation with proper headers."""
        snpeff_pipeline.generate_snippy_tab(
            self.vcf_file, self.out_tab, self.ann_file
        )
        self.assertTrue(os.path.exists(self.out_tab))
        with open(self.out_tab) as f:
            header = f.readline().strip()
        expected_header = "CHROM\tPOS\tTYPE\tREF\tALT\tFTYPE\tSTRAND\tNT_POS\tAA_POS\tEFFECT\tLOCUS_TAG\tGENE\tPRODUCT"
        self.assertEqual(header, expected_header)

    def test_generate_snippy_tab_annotates_upstream_and_cds(self):
        """Verify that upstream variants retain LOCUS_TAG, GENE, and PRODUCT annotations, with empty FTYPE."""
        gbk_content = (
            "LOCUS       chr1 1000 bp\n"
            "FEATURES             Location/Qualifiers\n"
            "     CDS             200..500\n"
            "                     /gene=\"testA\"\n"
            "                     /locus_tag=\"TAG_001\"\n"
            "                     /product=\"test protein\"\n"
            "//\n"
        )
        with open(self.ann_file, "w") as f:
            f.write(gbk_content)

        vcf_content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t50\t.\tC\tG\t.\t.\tANN=G|upstream_gene_variant|MODIFIER|testA|TAG_001|transcript|TAG_001|protein_coding||c.-150C>G|||||150|\n"
            "chr1\t250\t.\tA\tT\t.\t.\tANN=T|missense_variant|MODERATE|testA|TAG_001|transcript|TAG_001|protein_coding||c.51A>T|p.Lys17Asn|51/301|17/100|||\n"
        )
        with open(self.vcf_file, "w") as f:
            f.write(vcf_content)

        snpeff_pipeline.generate_snippy_tab(
            self.vcf_file, self.out_tab, self.ann_file
        )
        self.assertTrue(os.path.exists(self.out_tab))
        with open(self.out_tab) as f:
            lines = [l.strip().split("\t") for l in f if l.strip()]

        # Header + 2 variant rows
        self.assertEqual(len(lines), 3)

        # Row 1: upstream_gene_variant (FTYPE is empty, has TAG_001, testA, test protein)
        upstream_row = lines[1]
        self.assertEqual(upstream_row[0], "chr1")
        self.assertEqual(upstream_row[1], "50")
        self.assertEqual(upstream_row[5], "")  # FTYPE is empty for non-overlapping variant
        self.assertEqual(upstream_row[6], "+")  # STRAND from target feature
        self.assertIn("upstream_gene_variant", upstream_row[9])
        self.assertEqual(upstream_row[10], "TAG_001")
        self.assertEqual(upstream_row[11], "testA")
        self.assertEqual(upstream_row[12], "test protein")

        # Row 2: missense_variant (FTYPE is CDS)
        cds_row = lines[2]
        self.assertEqual(cds_row[0], "chr1")
        self.assertEqual(cds_row[1], "250")
        self.assertEqual(cds_row[5], "CDS")  # FTYPE is CDS for variant within CDS
        self.assertEqual(cds_row[6], "+")
        self.assertIn("missense_variant", cds_row[9])
        self.assertEqual(cds_row[10], "TAG_001")
        self.assertEqual(cds_row[11], "testA")
        self.assertEqual(cds_row[12], "test protein")



class TestHandleRun(unittest.TestCase):
    """Unit test for the 'run' handler logic without file I/O."""
    @patch("snpeff_pipeline.run_annotation")
    @patch("snpeff_pipeline.build_custom_db")
    def test_handle_run(self, mock_db_build, mock_annotate):
        args = snpeff_pipeline.parse_args([
            "run", "ref.fasta", "ref.gbk", "./input", "list.txt", "./output"
        ])
        build_config = DbBuildConfig.from_args(args)
        ann_config = AnnotationConfig.from_args(args)
        snpeff_pipeline.handle_run(build_config, ann_config)
        
        mock_db_build.assert_called_once_with(build_config)
        mock_annotate.assert_called_once_with(ann_config)


class TestEndToEndPipeline(unittest.TestCase):
    """End-to-End test for the 'run' subcommand pipeline with mock external calls."""
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.ref_fa = os.path.join(self.temp_dir, "ref_genome.fasta")
        with open(self.ref_fa, "w") as f:
            f.write(">ref_locus\nATGC\n")

        self.annotation_file = os.path.join(self.temp_dir, "ref_genome.gbk")
        with open(self.annotation_file, "w") as f:
            f.write("LOCUS       ref_locus 500 bp\n")

        self.input_dir = os.path.join(self.temp_dir, "vcfs")
        self.strain_dir = os.path.join(self.input_dir, "strain_01")
        os.makedirs(self.strain_dir)

        self.vcf_file = os.path.join(self.strain_dir, "snps.vcf")
        with open(self.vcf_file, "w") as f:
            f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nref_locus\t1\t.\tA\tT\t.\t.\t.\n")

        self.strain_list = os.path.join(self.temp_dir, "strain_list.txt")
        with open(self.strain_list, "w") as f:
            f.write("strain_01\n")

        self.output_dir = os.path.join(self.temp_dir, "results")
        self.db_dir = os.path.join(self.temp_dir, "built_db")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    @patch("snpeff_pipeline.get_base_config_path", return_value=None)
    @patch("snpeff_pipeline.generate_snippy_tab")
    @patch("snpeff_pipeline.run_command")
    def test_run_e2e_pipeline(self, mock_run_command, mock_snippy_tab, mock_get_base_config):
        """Test full pipeline execution flow (DB Build -> Strain Annotation)."""
        args = snpeff_pipeline.parse_args([
            "run", 
            self.ref_fa,
            self.annotation_file, 
            self.input_dir, 
            self.strain_list, 
            self.output_dir, 
            "--db_dir", self.db_dir
        ])
        build_config = DbBuildConfig.from_args(args)
        ann_config = AnnotationConfig.from_args(args)

        # Execute full pipeline (handle_run)
        snpeff_pipeline.handle_run(build_config, ann_config)
        # 1 call for snpEff build, 1 call for snpEff ann
        self.assertEqual(mock_run_command.call_count, 2)
        mock_snippy_tab.assert_called_once()


class TestE2ESubprocessPipeline(unittest.TestCase):
    """
    True End-to-End tests invoking snpeff_pipeline.py via external subprocess.
    Tests run with real test data and actual bioinformatics tools (snpEff, bcftools, snippy-vcf_to_tab).
    Output files are saved in timestamped date folders under test_output/YYYYMMDD_HHMMSS/
    to allow comparison between previous and current test runs.
    """
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    PIPELINE_SCRIPT = os.path.abspath(snpeff_pipeline.__file__)
    TEST_DATA_DIR = os.path.join(BASE_DIR, "test_data")
    TEST_OUTPUT_BASE = os.path.join(BASE_DIR, "test_output")

    @classmethod
    def setUpClass(cls):
        # Create timestamped date folder for this test suite execution
        cls.run_timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        cls.run_dir = os.path.join(cls.TEST_OUTPUT_BASE, cls.run_timestamp)
        os.makedirs(cls.run_dir, exist_ok=True)

        cls.ref_fasta = os.path.join(cls.TEST_DATA_DIR, "ref.fasta")
        cls.ref_auto = os.path.join(cls.TEST_DATA_DIR, "ref_auto.fasta")
        cls.ref_explicit = os.path.join(cls.TEST_DATA_DIR, "ref_explicit.fasta")
        cls.ref_gbk = os.path.join(cls.TEST_DATA_DIR, "ref.gbk")
        cls.ref_multi_gbk = os.path.join(cls.TEST_DATA_DIR, "ref_multi_contig.gbk")
        cls.ref_diff_locus_version = os.path.join(cls.TEST_DATA_DIR, "ref_diff_locus_version.gbk")
        cls.ref_gff = os.path.join(cls.TEST_DATA_DIR, "ref.gff")
        cls.ref_gtf = os.path.join(cls.TEST_DATA_DIR, "ref.gtf")
        cls.input_vcf = os.path.join(cls.TEST_DATA_DIR, "input_vcf")
        cls.input_vcf_auto = os.path.join(cls.TEST_DATA_DIR, "input_vcf_auto_map")
        cls.input_vcf_explicit = os.path.join(cls.TEST_DATA_DIR, "input_vcf_explicit_map")
        cls.chr_map_explicit = os.path.join(cls.TEST_DATA_DIR, "chr_map_explicit.txt")
        cls.strain_single = os.path.join(cls.TEST_DATA_DIR, "strain_list_single.txt")
        cls.strain_multi = os.path.join(cls.TEST_DATA_DIR, "strain_list_multi.txt")
        cls.strain_missing = os.path.join(cls.TEST_DATA_DIR, "strain_list_missing.txt")

        cls.expected = os.path.join(cls.TEST_DATA_DIR, "expected")

    def run_pipeline(self, *args) -> subprocess.CompletedProcess:
        """Helper to invoke snpeff_pipeline.py via subprocess."""
        cmd = [sys.executable, self.PIPELINE_SCRIPT] + list(args)
        return subprocess.run(cmd, capture_output=True, text=True)

    def assert_files_match(
        self,
        actual_file: str,
        expected_file: str,
        chr_remap: Optional[Dict[str, str]] = None
    ):
        """
        Compares actual and expected files line-by-line.
        Excludes:
        - Timestamp lines: 'Date ,' in summary CSV
        - '##SnpEffCmd=' in .ann.vcf
        - '##reference=' in .ann.vcf (when comparing cross-reference rename tests)
        - 'Command_line_arguments ,' in summary CSV
        All other lines including headers and variant data are strictly compared.
        """
        self.assertTrue(os.path.exists(actual_file), f"Actual file missing: {actual_file}")
        self.assertTrue(os.path.exists(expected_file), f"Expected file missing: {expected_file}")

        with open(actual_file, "r") as f_act, open(expected_file, "r") as f_exp:
            act_lines = [
                line for line in f_act
                if not line.startswith("##SnpEffCmd")
                and not line.startswith("Command_line_arguments")
                and not line.startswith("Date ,")
                and not line.startswith("##reference=")
            ]
            exp_lines = [
                line for line in f_exp
                if not line.startswith("##SnpEffCmd")
                and not line.startswith("Command_line_arguments")
                and not line.startswith("Date ,")
                and not line.startswith("##reference=")
            ]

        if chr_remap:
            new_act_lines = []
            for line in act_lines:
                for k, v in chr_remap.items():
                    line = line.replace(k, v)
                new_act_lines.append(line)
            act_lines = new_act_lines

        self.assertEqual(
            len(act_lines), len(exp_lines),
            f"Line count mismatch between {actual_file} ({len(act_lines)}) and {expected_file} ({len(exp_lines)})"
        )
        for i, (act, exp) in enumerate(zip(act_lines, exp_lines), 1):
            self.assertEqual(act, exp, f"Mismatch at line {i} in {actual_file} vs {expected_file}")

    def assert_output_directory_matches(
        self,
        actual_dir: str,
        expected_dir: str,
        strain: str = "strain_A",
        chr_remap: Optional[Dict[str, str]] = None
    ):
        """Verify all core result files match the expected directory."""
        files_to_check = [
            f"{strain}.ann.vcf",
            f"{strain}.snps.tab",
            f"{strain}_summary.csv",
            f"{strain}_summary.genes.txt"
        ]
        for fname in files_to_check:
            self.assert_files_match(
                os.path.join(actual_dir, fname),
                os.path.join(expected_dir, fname),
                chr_remap=chr_remap
            )

    def test_01_run_pipeline_gbk_matches_expected(self):
        """Test 1: 'run' subcommand with GenBank annotation matches expected output."""
        case_dir = os.path.join(self.run_dir, "test_01_run_gbk")
        out_dir = os.path.join(case_dir, "output")
        db_dir = os.path.join(case_dir, "db")

        proc = self.run_pipeline(
            "run", self.ref_fasta, self.ref_gbk, self.input_vcf, self.strain_single,
            out_dir, "--db_dir", db_dir
        )
        self.assertEqual(proc.returncode, 0, f"Pipeline failed: {proc.stderr}")
        # Matching chromosome names -> no auto-mapping needed
        self.assertNotIn("Auto-mapped single chromosome", proc.stderr + proc.stdout)
        self.assert_output_directory_matches(out_dir, self.expected)

    def test_02_run_pipeline_gff_matches_expected(self):
        """Test 2: 'run' subcommand with GFF3 annotation matches expected output."""
        case_dir = os.path.join(self.run_dir, "test_02_run_gff")
        out_dir = os.path.join(case_dir, "output")
        db_dir = os.path.join(case_dir, "db")

        proc = self.run_pipeline(
            "run", self.ref_fasta, self.ref_gff, self.input_vcf, self.strain_single,
            out_dir, "--db_dir", db_dir
        )
        self.assertEqual(proc.returncode, 0, f"Pipeline failed: {proc.stderr}")
        self.assert_output_directory_matches(out_dir, self.expected)

    def test_03_run_pipeline_gtf_matches_expected(self):
        """Test 3: 'run' subcommand with GTF annotation matches expected output."""
        case_dir = os.path.join(self.run_dir, "test_03_run_gtf")
        out_dir = os.path.join(case_dir, "output")
        db_dir = os.path.join(case_dir, "db")

        proc = self.run_pipeline(
            "run", self.ref_fasta, self.ref_gtf, self.input_vcf, self.strain_single,
            out_dir, "--db_dir", db_dir
        )
        self.assertEqual(proc.returncode, 0, f"Pipeline failed: {proc.stderr}")
        self.assert_output_directory_matches(out_dir, self.expected)

    def test_04_db_build_then_annotate_matches_expected(self):
        """Test 4: Separate 'db-build' then 'annotate' execution matches expected output."""
        case_dir = os.path.join(self.run_dir, "test_04_db_build_then_annotate")
        out_dir = os.path.join(case_dir, "output")
        db_dir = os.path.join(case_dir, "db")

        # Step 1: db-build
        proc_build = self.run_pipeline(
            "db-build", self.ref_fasta, self.ref_gbk, "--db_dir", db_dir
        )
        self.assertEqual(proc_build.returncode, 0, f"db-build failed: {proc_build.stderr}")

        # Step 2: annotate
        proc_ann = self.run_pipeline(
            "annotate", self.input_vcf, self.strain_single, out_dir, "--db_dir", db_dir
        )
        self.assertEqual(proc_ann.returncode, 0, f"annotate failed: {proc_ann.stderr}")
        self.assert_output_directory_matches(out_dir, self.expected)

    def test_05_auto_chromosome_rename(self):
        """Test 5: Auto chromosome adaptation from differing annotation CHROM (BA000007.3) to reference FASTA CHROM (Sakai_BA000007)."""
        case_dir = os.path.join(self.run_dir, "test_05_auto_chr_rename")
        out_dir = os.path.join(case_dir, "output")
        db_dir = os.path.join(case_dir, "db")

        # Use ref_auto (header: Sakai_BA000007) and input_vcf_auto (CHROM: Sakai_BA000007), annotation is ref_gbk (BA000007.3)
        proc = self.run_pipeline(
            "run", self.ref_auto, self.ref_gbk, self.input_vcf_auto, self.strain_single,
            out_dir, "--db_dir", db_dir
        )
        self.assertEqual(proc.returncode, 0, f"Pipeline failed: {proc.stderr}")
        self.assertIn("Auto-mapped single chromosome for DB build: 'BA000007.3' -> 'Sakai_BA000007'", proc.stderr + proc.stdout)
        self.assert_output_directory_matches(
            out_dir, self.expected, chr_remap={"Sakai_BA000007": "BA000007.3"}
        )

    def test_06_explicit_chr_map(self):
        """Test 6: Explicit --chr_map adapts annotation CHROM (BA000007.3) to reference FASTA CHROM (chr_explicit)."""
        case_dir = os.path.join(self.run_dir, "test_06_explicit_chr_map")
        out_dir = os.path.join(case_dir, "output")
        db_dir = os.path.join(case_dir, "db")

        # Use ref_explicit (header: chr_explicit), input_vcf_explicit (CHROM: chr_explicit), mapped via chr_map_explicit (BA000007.3 -> chr_explicit)
        proc = self.run_pipeline(
            "run", self.ref_explicit, self.ref_gbk, self.input_vcf_explicit, self.strain_single,
            out_dir, "--db_dir", db_dir, "--chr_map", self.chr_map_explicit
        )
        self.assertEqual(proc.returncode, 0, f"Pipeline failed: {proc.stderr}")
        self.assertIn("Using explicit chromosome mapping for DB build", proc.stderr + proc.stdout)
        self.assert_output_directory_matches(
            out_dir, self.expected, chr_remap={"chr_explicit": "BA000007.3", "explicit": "BA000007.3"}
        )

    def test_07_multiple_strains_batch(self):
        """Test 7: Batch processing multiple strains (strain_A + strain_B)."""
        case_dir = os.path.join(self.run_dir, "test_07_multiple_strains")
        out_dir = os.path.join(case_dir, "output")
        db_dir = os.path.join(case_dir, "db")

        proc = self.run_pipeline(
            "run", self.ref_fasta, self.ref_gbk, self.input_vcf, self.strain_multi,
            out_dir, "--db_dir", db_dir
        )
        self.assertEqual(proc.returncode, 0, f"Pipeline failed: {proc.stderr}")
        self.assert_output_directory_matches(out_dir, self.expected, strain="strain_A")
        self.assert_output_directory_matches(out_dir, self.expected, strain="strain_B")

    def test_08_skip_missing_strain(self):
        """Test 8: Missing strain VCF is gracefully skipped without breaking the pipeline."""
        case_dir = os.path.join(self.run_dir, "test_08_skip_missing_strain")
        out_dir = os.path.join(case_dir, "output")
        db_dir = os.path.join(case_dir, "db")

        proc = self.run_pipeline(
            "run", self.ref_fasta, self.ref_gbk, self.input_vcf, self.strain_missing,
            out_dir, "--db_dir", db_dir
        )
        self.assertEqual(proc.returncode, 0, f"Pipeline failed: {proc.stderr}")
        self.assert_output_directory_matches(out_dir, self.expected, strain="strain_A")
        self.assertFalse(os.path.exists(os.path.join(out_dir, "nonexistent_strain.ann.vcf")))
        self.assertIn("Input VCF not found for nonexistent_strain", proc.stderr + proc.stdout)

    def test_09_annotate_without_db_fails(self):
        """Test 9: 'annotate' command fails when pointing to a non-existent database directory."""
        case_dir = os.path.join(self.run_dir, "test_09_annotate_without_db")
        out_dir = os.path.join(case_dir, "output")
        nonexistent_db = os.path.join(case_dir, "no_such_db")

        proc = self.run_pipeline(
            "annotate", self.input_vcf, self.strain_single, out_dir, "--db_dir", nonexistent_db
        )
        self.assertGreaterEqual(proc.returncode, 1, f"Expected non-zero exit code (>= 1), got {proc.returncode}")
        self.assertEqual(proc.returncode, 1)
        self.assertIn("Database directory not found", proc.stderr + proc.stdout)

    def test_10_reverse_chr_map_fails(self):
        """Test 10: Explicit --chr_map in reverse order ([ref_fasta] [ann]) fails with an error."""
        case_dir = os.path.join(self.run_dir, "test_10_reverse_chr_map")
        out_dir = os.path.join(case_dir, "output")
        db_dir = os.path.join(case_dir, "db")
        os.makedirs(case_dir, exist_ok=True)

        reverse_map_file = os.path.join(case_dir, "reverse_map.txt")
        with open(reverse_map_file, "w") as f:
            # ref_fasta CHROM first, annotation CHROM second (reversed order!)
            f.write("chr_explicit\tBA000007.3\n")

        proc = self.run_pipeline(
            "run", self.ref_explicit, self.ref_gbk, self.input_vcf_explicit, self.strain_single,
            out_dir, "--db_dir", db_dir, "--chr_map", reverse_map_file
        )
        self.assertGreaterEqual(proc.returncode, 1, f"Expected non-zero exit code, got {proc.returncode}")
        self.assertIn("not found in annotation file chromosomes", proc.stderr + proc.stdout)

    def test_11_gbk_diff_locus_version(self):
        """Test 11: Real GenBank file where LOCUS ('BA000007') differs from VERSION ('BA000007.3')."""
        case_dir = os.path.join(self.run_dir, "test_11_gbk_diff_locus_version")
        out_dir = os.path.join(case_dir, "output")
        db_dir = os.path.join(case_dir, "db")

        # ref_auto has CHROM Sakai_BA000007, ref_diff_locus_version has LOCUS BA000007 and VERSION BA000007.3
        proc = self.run_pipeline(
            "run", self.ref_auto, self.ref_diff_locus_version, self.input_vcf_auto, self.strain_single,
            out_dir, "--db_dir", db_dir
        )
        self.assertEqual(proc.returncode, 0, f"Pipeline failed: {proc.stderr}")
        self.assertIn("Auto-mapped single chromosome for DB build: 'BA000007.3' -> 'Sakai_BA000007'", proc.stderr + proc.stdout)
        self.assert_output_directory_matches(
            out_dir, self.expected, strain="strain_A", chr_remap={"Sakai_BA000007": "BA000007.3"}
        )

    def test_12_partial_chr_map_multi_contig(self):
        """Test 12: Multi-contig annotation file with only the matching chromosome specified in --chr_map."""
        case_dir = os.path.join(self.run_dir, "test_12_partial_chr_map")
        out_dir = os.path.join(case_dir, "output")
        db_dir = os.path.join(case_dir, "db")

        # ref_multi_gbk has ['BA000007.3', 'BA000008.1'], chr_map_explicit only maps BA000007.3 -> chr_explicit
        proc = self.run_pipeline(
            "run", self.ref_explicit, self.ref_multi_gbk, self.input_vcf_explicit, self.strain_single,
            out_dir, "--db_dir", db_dir, "--chr_map", self.chr_map_explicit
        )
        self.assertEqual(proc.returncode, 0, f"Pipeline failed: {proc.stderr}")
        self.assertIn("Using explicit chromosome mapping for DB build", proc.stderr + proc.stdout)
        self.assert_output_directory_matches(
            out_dir, self.expected, strain="strain_A", chr_remap={"chr_explicit": "BA000007.3", "explicit": "BA000007.3"}
        )

    def test_13_chr_map_identical_names_allowed(self):
        """Test 13: --chr_map with identical source and target names is allowed."""
        case_dir = os.path.join(self.run_dir, "test_13_identical_names")
        out_dir = os.path.join(case_dir, "output")
        db_dir = os.path.join(case_dir, "db")
        os.makedirs(case_dir, exist_ok=True)

        map_file = os.path.join(case_dir, "map_identical.txt")
        with open(map_file, "w") as f:
            f.write("BA000007.3\tBA000007.3\n")

        proc = self.run_pipeline(
            "run", self.ref_fasta, self.ref_gbk, self.input_vcf, self.strain_single,
            out_dir, "--db_dir", db_dir, "--chr_map", map_file
        )
        self.assertEqual(proc.returncode, 0, f"Pipeline failed: {proc.stderr}")
        self.assert_output_directory_matches(out_dir, self.expected)

    def test_14_chr_map_duplicate_source_fails(self):
        """Test 14: --chr_map with duplicate annotation chromosome specified fails."""
        case_dir = os.path.join(self.run_dir, "test_14_duplicate_source")
        out_dir = os.path.join(case_dir, "output")
        db_dir = os.path.join(case_dir, "db")
        os.makedirs(case_dir, exist_ok=True)

        map_file = os.path.join(case_dir, "map_dup_source.txt")
        with open(map_file, "w") as f:
            f.write("BA000007.3\tchr_explicit\nBA000007.3\tchr_other\n")

        proc = self.run_pipeline(
            "run", self.ref_explicit, self.ref_gbk, self.input_vcf_explicit, self.strain_single,
            out_dir, "--db_dir", db_dir, "--chr_map", map_file
        )
        self.assertGreaterEqual(proc.returncode, 1)
        self.assertIn("is specified more than once", proc.stderr + proc.stdout)

    def test_15_chr_map_duplicate_target_fails(self):
        """Test 15: --chr_map mapping multiple annotation chromosomes to the same reference chromosome fails."""
        case_dir = os.path.join(self.run_dir, "test_15_duplicate_target")
        out_dir = os.path.join(case_dir, "output")
        db_dir = os.path.join(case_dir, "db")
        os.makedirs(case_dir, exist_ok=True)

        map_file = os.path.join(case_dir, "map_dup_target.txt")
        with open(map_file, "w") as f:
            f.write("BA000007.3\tchr_explicit\nBA000008.1\tchr_explicit\n")

        proc = self.run_pipeline(
            "run", self.ref_explicit, self.ref_multi_gbk, self.input_vcf_explicit, self.strain_single,
            out_dir, "--db_dir", db_dir, "--chr_map", map_file
        )
        self.assertGreaterEqual(proc.returncode, 1)
        self.assertIn("is already mapped from another chromosome", proc.stderr + proc.stdout)

    def test_16_chr_map_target_not_in_ref_fails(self):
        """Test 16: --chr_map with target chromosome not in reference FASTA fails."""
        case_dir = os.path.join(self.run_dir, "test_16_target_not_in_ref")
        out_dir = os.path.join(case_dir, "output")
        db_dir = os.path.join(case_dir, "db")
        os.makedirs(case_dir, exist_ok=True)

        map_file = os.path.join(case_dir, "map_not_in_ref.txt")
        with open(map_file, "w") as f:
            f.write("BA000007.3\tunknown_ref_name\n")

        proc = self.run_pipeline(
            "run", self.ref_explicit, self.ref_gbk, self.input_vcf_explicit, self.strain_single,
            out_dir, "--db_dir", db_dir, "--chr_map", map_file
        )
        self.assertGreaterEqual(proc.returncode, 1)
        self.assertIn("not found in reference FASTA chromosomes", proc.stderr + proc.stdout)


if __name__ == "__main__":
    unittest.main()



