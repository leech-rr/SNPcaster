import os
import sys
import pytest
import pandas as pd
from unittest.mock import patch, MagicMock
from merge_amrfinder import (
    get_sort_key, get_group_name, merge_amrfinder_files,
    load_amrfinder_files, generate_summary_tsv, generate_merged_tsv
)

# --- Unit Tests for get_sort_key ---

@pytest.mark.parametrize(
    "scope, type_str, expected",
    [
        ("core", "AMR", (0, "core", 0, "amr")),
        ("plus", "AMR", (1, "plus", 0, "amr")),
        ("plus", "VIRULENCE", (1, "plus", 1, "virulence")),
        ("plus", "STRESS", (1, "plus", 2, "stress")),
        ("core", "STRESS", (0, "core", 2, "stress")),
        ("other", "AMR", (2, "other", 0, "amr")),
        ("plus", "other", (1, "plus", 3, "other")),
    ]
)
def test_get_sort_key(scope, type_str, expected):
    """Test get_sort_key with various combinations of Scope and Type."""
    assert get_sort_key(scope, type_str) == expected


# --- Unit Tests for get_group_name ---

@pytest.mark.parametrize(
    "scope, type_str, expected",
    [
        ("core", "AMR", "core_AMR"),
        ("plus", "AMR", "plus_AMR"),
        ("plus", "VIRULENCE", "plus_virulence"),
        ("plus", "STRESS", "plus_stress"),
        ("core", "STRESS", "core_stress"),
        ("other", "AMR", "other_AMR"),
        ("plus", "unknown", "plus_unknown"),
    ]
)
def test_get_group_name(scope, type_str, expected):
    """Test get_group_name formatting."""
    assert get_group_name(scope, type_str) == expected


# --- Fixtures for Mock Data ---

@pytest.fixture
def mock_normal_data(tmp_path):
    """
    Fixture to set up valid strain list and valid TSV files.
    All files are created inside the pytest managed `tmp_path` directory,
    which is automatically cleaned up by pytest after test execution.
    """
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    
    # 1. Create list_amrfinder.txt
    list_path = tmp_path / "list_amrfinder.txt"
    strains = ["StrainA", "StrainB"]
    list_path.write_text("\n".join(strains), encoding="utf-8")
    
    # 2. Create TSV for StrainA
    tsv_a = input_dir / "StrainA_plus_out.tsv"
    df_a = pd.DataFrame({
        "Protein id": ["p1", "p2"],
        "Element symbol": ["geneX", "geneY"],
        "Scope": ["core", "plus"],
        "Type": ["AMR", "VIRULENCE"],
        "Extra Col": ["val1", "val2"]
    })
    df_a.to_csv(tsv_a, sep="\t", index=False)
    
    # 3. Create TSV for StrainB
    tsv_b = input_dir / "StrainB_plus_out.tsv"
    df_b = pd.DataFrame({
        "Protein id": ["p3"],
        "Element symbol": ["geneX"],
        "Scope": ["core"],
        "Type": ["AMR"]
    })
    df_b.to_csv(tsv_b, sep="\t", index=False)
    
    return {
        "input_dir": str(input_dir),
        "list_path": str(list_path),
        "output_dir": str(output_dir),
        "strains": strains
    }


# --- Integration Tests for merge_amrfinder_files ---

def test_merge_amrfinder_files_normal(mock_normal_data):
    """Test normal end-to-end execution of merge_amrfinder_files."""
    input_dir = mock_normal_data["input_dir"]
    list_path = mock_normal_data["list_path"]
    output_dir = mock_normal_data["output_dir"]
    
    # Run the merger
    merge_amrfinder_files(input_dir, list_path, output_dir)
    
    # Assert output file exists
    expected_output_path = os.path.join(output_dir, "list_amrfinder_summary.tsv")
    assert os.path.exists(expected_output_path)
    
    # Read output and verify structure
    df_res = pd.read_csv(expected_output_path, sep="\t")
    
    # Expected columns: Strain, core_AMR, geneX (sorted), plus_virulence, geneY
    expected_columns = ["Strain", "core_AMR", "geneX", "plus_virulence", "geneY"]
    assert list(df_res.columns) == expected_columns
    
    # Expected rows
    assert len(df_res) == 2
    
    row_a = df_res[df_res["Strain"] == "StrainA"].iloc[0]
    assert row_a["core_AMR"] == "|"
    assert row_a["geneX"] == 1
    assert row_a["plus_virulence"] == "|"
    assert row_a["geneY"] == 1
    
    row_b = df_res[df_res["Strain"] == "StrainB"].iloc[0]
    assert row_b["core_AMR"] == "|"
    assert row_b["geneX"] == 1
    assert row_b["plus_virulence"] == "|"
    assert row_b["geneY"] == 0

    # Assert output file exists for merged concatenation
    expected_merged_path = os.path.join(output_dir, "list_amrfinder_merged.tsv")
    assert os.path.exists(expected_merged_path)
    
    # Read merged output and verify structure and order
    df_merged = pd.read_csv(expected_merged_path, sep="\t")
    
    # Expected columns in merged file: Strain at start, then the original columns
    assert df_merged.columns[0] == "Strain"
    assert "Protein id" in df_merged.columns
    assert "Element symbol" in df_merged.columns
    assert "Scope" in df_merged.columns
    assert "Type" in df_merged.columns
    
    # Strain order and row order verification:
    # StrainA has 2 rows (p1, p2), StrainB has 1 row (p3)
    # Total row count should be 3
    assert len(df_merged) == 3
    
    # Row 0: StrainA, p1, geneX, core, AMR, val1
    assert df_merged.iloc[0]["Strain"] == "StrainA"
    assert df_merged.iloc[0]["Protein id"] == "p1"
    assert df_merged.iloc[0]["Element symbol"] == "geneX"
    assert df_merged.iloc[0]["Scope"] == "core"
    assert df_merged.iloc[0]["Type"] == "AMR"
    assert df_merged.iloc[0]["Extra Col"] == "val1"
    
    # Row 1: StrainA, p2, geneY, plus, VIRULENCE, val2
    assert df_merged.iloc[1]["Strain"] == "StrainA"
    assert df_merged.iloc[1]["Protein id"] == "p2"
    assert df_merged.iloc[1]["Element symbol"] == "geneY"
    assert df_merged.iloc[1]["Scope"] == "plus"
    assert df_merged.iloc[1]["Type"] == "VIRULENCE"
    assert df_merged.iloc[1]["Extra Col"] == "val2"
    
    # Row 2: StrainB, p3, geneX, core, AMR, NaN (missing column Extra Col in StrainB)
    assert df_merged.iloc[2]["Strain"] == "StrainB"
    assert df_merged.iloc[2]["Protein id"] == "p3"
    assert df_merged.iloc[2]["Element symbol"] == "geneX"
    assert df_merged.iloc[2]["Scope"] == "core"
    assert df_merged.iloc[2]["Type"] == "AMR"
    # Column "Extra Col" was only in StrainA, so pd.concat should yield NaN for StrainB
    assert pd.isna(df_merged.iloc[2]["Extra Col"])


# --- Integration Tests for Error Scenarios ---

@pytest.mark.parametrize(
    "error_type, expected_exception, error_msg_substring",
    [
        ("missing_tsv", FileNotFoundError, "not found"),
        ("empty_tsv", ValueError, "is empty"),
        ("missing_columns", KeyError, "Missing required columns"),
        ("empty_list", ValueError, "No strain IDs found"),
        ("missing_list", FileNotFoundError, "does not exist"),
    ]
)
def test_merge_errors(tmp_path, error_type, expected_exception, error_msg_substring):
    """
    Test merge_amrfinder_files error scenarios using parameterized inputs.
    All files are created inside the pytest managed `tmp_path` directory,
    which is automatically cleaned up by pytest after test execution.
    """
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    
    list_path = tmp_path / "list_amrfinder.txt"
    strains = ["StrainA"]
    
    if error_type == "missing_list":
        # Do not create list file
        pass
    elif error_type == "empty_list":
        # Create empty list file
        list_path.write_text("", encoding="utf-8")
    else:
        # Create valid list file
        list_path.write_text("\n".join(strains), encoding="utf-8")
        
        tsv_path = input_dir / "StrainA_plus_out.tsv"
        if error_type == "missing_tsv":
            # Do not create TSV file
            pass
        elif error_type == "empty_tsv":
            # Create empty file
            tsv_path.write_text("", encoding="utf-8")
        elif error_type == "missing_columns":
            # Create TSV with missing column 'Type'
            df = pd.DataFrame({
                "Element symbol": ["geneX"],
                "Scope": ["core"]
            })
            df.to_csv(tsv_path, sep="\t", index=False)
            
    # Execute and verify expected exceptions
    with pytest.raises(expected_exception) as excinfo:
        merge_amrfinder_files(str(input_dir), str(list_path), str(output_dir))
        
    assert error_msg_substring in str(excinfo.value)


# --- Unit Tests for main() ---

def test_main_success():
    """Test main() parses CLI arguments correctly and calls merge_amrfinder_files."""
    test_args = ["merge_amrfinder.py", "mock_input_dir", "mock_list.txt", "mock_output_dir"]
    with patch("sys.argv", test_args), \
         patch("merge_amrfinder.merge_amrfinder_files") as mock_merge:
        
        from merge_amrfinder import main
        main()
        
        # Verify merge_amrfinder_files was called with the correct args
        mock_merge.assert_called_once_with("mock_input_dir", "mock_list.txt", "mock_output_dir")


@pytest.mark.parametrize(
    "exception_to_raise",
    [
        FileNotFoundError("mocked file not found error"),
        ValueError("mocked value error"),
        OSError("mocked OS error"),
        Exception("mocked unexpected error"),
    ]
)
def test_main_error_handling(exception_to_raise):
    """Test main() catches exceptions and exits with exit code 1."""
    test_args = ["merge_amrfinder.py", "mock_input_dir", "mock_list.txt", "mock_output_dir"]
    
    with patch("sys.argv", test_args), \
         patch("merge_amrfinder.merge_amrfinder_files", side_effect=exception_to_raise), \
         patch("sys.stderr", new_callable=MagicMock) as mock_stderr:
        
        from merge_amrfinder import main
        
        # main() should call sys.exit(), raising SystemExit
        with pytest.raises(SystemExit) as excinfo:
            main()
            
        # Verify exit code is 1
        assert excinfo.value.code == 1
        
        # Verify that the error message was printed to stderr
        mock_stderr.write.assert_called()
        printed_msg = "".join(call.args[0] for call in mock_stderr.write.call_args_list)
        assert str(exception_to_raise) in printed_msg


# --- Unit Tests for newly extracted functions ---

def test_load_amrfinder_files(mock_normal_data):
    """Test load_amrfinder_files to ensure it preserves the exact order of strains and loads raw lines."""
    input_dir = mock_normal_data["input_dir"]
    list_path = mock_normal_data["list_path"]
    
    list_name, strain_raw_pairs = load_amrfinder_files(input_dir, list_path)
    
    assert list_name == "list_amrfinder"
    
    # Verify exact insertion order in list of tuples
    assert len(strain_raw_pairs) == 2
    
    assert strain_raw_pairs[0][0] == "StrainA"
    assert isinstance(strain_raw_pairs[0][1], list)
    # StrainA has header + 2 data rows = 3 lines total
    assert len(strain_raw_pairs[0][1]) == 3
    
    assert strain_raw_pairs[1][0] == "StrainB"
    assert isinstance(strain_raw_pairs[1][1], list)
    # StrainB has header + 1 data row = 2 lines total
    assert len(strain_raw_pairs[1][1]) == 2


def test_generate_merged_tsv(mock_normal_data):
    """Test generate_merged_tsv generates output in the correct directory, preserves order, and adds Strain column."""
    input_dir = mock_normal_data["input_dir"]
    list_path = mock_normal_data["list_path"]
    output_dir = mock_normal_data["output_dir"]
    
    list_name, strain_raw_pairs = load_amrfinder_files(input_dir, list_path)
    generate_merged_tsv(strain_raw_pairs, list_name, output_dir)
    
    expected_merged_path = os.path.join(output_dir, "list_amrfinder_merged.tsv")
    assert os.path.exists(expected_merged_path)
    
    df_merged = pd.read_csv(expected_merged_path, sep="\t")
    
    # Total row count, columns, and ordering verification
    assert len(df_merged) == 3
    assert list(df_merged["Strain"]) == ["StrainA", "StrainA", "StrainB"]
    assert list(df_merged["Protein id"]) == ["p1", "p2", "p3"]
