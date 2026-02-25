import subprocess
import os
import h5py
from pathlib import Path


def test_population_cli_runs(tmp_path):
    config_file = Path(__file__).parent / "test.ini"
    output_file = tmp_path / "test.h5"

    cmd = [
        "gwforge_population",
        "--config-file",
        str(config_file),
        "--output-file",
        str(output_file),
    ]

    result = subprocess.run(cmd, capture_output=True)
    assert result.returncode == 0
    assert os.path.exists(output_file)

    with h5py.File(output_file, "r") as f:
        assert len(f.keys()) > 0
        assert "mass_1" in f.keys()
        assert "mass_2" in f.keys()
