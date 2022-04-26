"""Functions for parsing measurements from raw Maud inputs."""

import pandas as pd

from maud.data_model.measurement_set import (
    EnzymeKnockout,
    Experiment,
    MeasurementSet,
    MeasurementType,
    PhosphorylationKnockout,
)


def parse_measurement_set(
    raw_measurement_table: pd.DataFrame, raw_biological_config: dict
) -> MeasurementSet:
    """Parse a measurements dataframe.

    :param measurement_table: result of running pd.read_csv on suitable file
    :param raw_biological_config: result of running toml.load on suitable file
    """
    experiments = [Experiment(**e) for e in raw_biological_config["experiment"]]
    y = {
        mt: raw_measurement_table.loc[
            lambda df: df["measurement_type"] == mt.value
        ]
        for mt in MeasurementType
    }
    enzyme_knockouts = (
        [
            EnzymeKnockout(
                experiment_id=eko["experiment_id"], enzyme_id=eko["enzyme_id"]
            )
            for eko in raw_biological_config["enzyme_knockout"]
        ]
        if "enzyme_knockout" in raw_biological_config.keys()
        else None
    )
    phosphorylation_knockouts = (
        [
            PhosphorylationKnockout(
                experiment_id=pko["experiment_id"], enzyme_id=pko["enzyme_id"]
            )
            for pko in raw_biological_config["phos_knockout"]
        ]
        if "phos_knockout" in raw_biological_config.keys()
        else None
    )
    return MeasurementSet(
        yconc=y[MeasurementType.MIC],
        yflux=y[MeasurementType.FLUX],
        yenz=y[MeasurementType.ENZYME],
        enzyme_knockouts=enzyme_knockouts,
        phosphorylation_knockouts=phosphorylation_knockouts,
        experiments=experiments,
    )
