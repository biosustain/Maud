"""Test functions in the cli module."""

from click.testing import CliRunner

from maud.cli import cli


def test_maud_help():
    """Test that running `maud --help` does not raise an error."""
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
