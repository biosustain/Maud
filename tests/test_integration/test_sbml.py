from pathlib import Path

from maud.sbml_compat import MaudModelBuilder, SBMLModelParser


data_dir = Path(__file__).parent.parent / "data"


def test_sbml_parser():
    parser = SBMLModelParser.from_file(data_dir / "BIOMD0000000639_urn.xml")
    model = parser.parse(MaudModelBuilder())
    assert model.model_id == "MODEL1602280001"
    assert len(model.metabolites) == 17
    assert len(model.reactions) == 10
