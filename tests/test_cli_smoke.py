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
    "gwforge_population_fisher",
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


def test_population_fisher_compare_writes_a_figure(tmp_path):
    """``--compare`` overlays saved forecasts without recomputing them.

    Exercised end to end because it is the path that produces the
    network-versus-network figure, and it has no other test: the plotting is
    covered in ``test_plotting.py`` and the Fisher in
    ``test_population_fisher.py``, but nothing else runs the two together
    through the CLI.
    """
    import numpy

    names = numpy.array(["gamma", "kappa"], dtype=object)
    paths = []
    for index, scale in enumerate((1.0, 4.0)):
        path = tmp_path / "forecast_{}.npz".format(index)
        numpy.savez(
            path,
            fisher=numpy.eye(2) / scale**2,
            covariance=numpy.eye(2) * scale**2,
            parameter_names=names,
            fiducial=numpy.array([2.7, 5.6]),
            snr_threshold=10.0,
            n_events=100,
            n_total=1000,
            condition_number=1.0,
        )
        paths.append(path)

    output = tmp_path / "comparison.pdf"
    result = subprocess.run(
        [
            sys.executable,
            str(BIN / "gwforge_population_fisher"),
            "--compare",
            "A={}".format(paths[0]),
            "B={}".format(paths[1]),
            "--output-file",
            str(output),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert output.exists()
    # The ratio column is what the figure is for; B is four times wider.
    assert "4" in result.stdout
