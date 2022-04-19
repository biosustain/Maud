"""Unit tests for sampling functions."""

import shutil
import tempfile
from pathlib import Path

from maud.cli import generate_predictions, simulate
from maud.plotting.plot_oos import plot_oos
from maud.plotting.plot_samples import plot_posteriors


home_dir = Path.cwd()
data_wd = home_dir / "tests" / "data"
n_samples = 5


def test_plot_samples():
    """Smoke test of plot_samples script.

    Ensures that the plot_samples script is functional
    and tests this by parsing a Maud output into the plotting
    function. Ensures name changes are updated in this script
    """
    linear_model_path = data_wd / "linear"
    tmp_dir = tempfile.mkdtemp(dir=data_wd)
    output_path = Path(tmp_dir)
    maud_output_path = Path(
        simulate(linear_model_path, output_path, n_samples)
    )
    plot_posteriors(maud_output_path, output_path)
    shutil.rmtree(tmp_dir)
    assert True


def test_plot_oos():
    """Smoke test of plot_oos script.

    Ensures that the plot_oos script is functional
    and tests this by parsing a Maud output into the plotting
    function. Ensures name changes are updated in this script
    """
    linear_model_path = data_wd / "linear"
    tmp_dir = tempfile.mkdtemp(dir=data_wd)
    output_path = Path(tmp_dir)
    maud_output_path = Path(
        simulate(linear_model_path, output_path, n_samples)
    )
    maud_oos_path = Path(
        generate_predictions(maud_output_path, linear_model_path, output_path)
    )
    plot_oos(maud_oos_path, output_path)
    shutil.rmtree(tmp_dir)
    assert True
