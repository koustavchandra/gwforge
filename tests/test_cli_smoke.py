"""Smoke tests: every CLI must at least import and print --help without crashing.

These would have caught the `split_odd_even` ImportError in gwforge_workflow, the
hardcoded interpreter shebangs, and the missing `import os` in gwforge_inject.
"""

import subprocess
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parent.parent / "bin"
SCRIPTS = [
    "gwforge_population",
    "gwforge_noise",
    "gwforge_inject",
    "gwforge_workflow",
    "gwforge_optimal_snr",
    "gwforge_estimate_psd",
    "gwforge_fisher",
]


@pytest.mark.parametrize("script", SCRIPTS)
def test_cli_help(script):
    # Invoke through the current interpreter so the test is independent of the
    # (non-portable) shebang and of whether console scripts are on PATH.
    result = subprocess.run(
        [sys.executable, str(BIN / script), "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"{script} --help failed:\n{result.stderr}"
    assert "usage" in result.stdout.lower()
