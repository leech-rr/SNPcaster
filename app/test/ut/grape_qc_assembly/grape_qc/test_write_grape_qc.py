import os
import tempfile
from unittest.mock import patch, MagicMock
from src.grape_qc_assembly.write_grape_qc import calculate_read_stats

def test_calculate_read_stats():
    with tempfile.TemporaryDirectory() as temp_dir:
        strain = "strain1"
        fastq1 = os.path.join(temp_dir, f"{strain}_1.fastq.gz")
        with open(fastq1, 'w') as f:
            f.write("dummy")
            
        with patch('src.grape_qc_assembly.write_grape_qc.subprocess.run') as mock_run:
            mock_result = MagicMock()
            mock_result.stdout = "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\nstrain1_1.fastq.gz\tFASTQ\tDNA\t100\t15000\t150\t150\t150\n"
            mock_run.return_value = mock_result
            
            reads, total_len = calculate_read_stats(temp_dir, [f"{strain}_1.fastq.gz"])
            
            assert reads == 100
            assert total_len == 15000
            print("test_calculate_read_stats passed")

if __name__ == "__main__":
    test_calculate_read_stats()
    print("All tests passed!")
