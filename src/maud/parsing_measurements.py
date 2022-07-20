"""Functions for parsing measurements from raw Maud inputs."""

import pandas as pd

from maud.data_model.measurement_set import (
    EnzymeKnockout,
    Experiment,
    MeasurementSet,
    MeasurementType,
    PhosphorylationModifyingEnzymeKnockout,
)


def parse_measurement_set(
    raw_measurement_table: pd.DataFrame, raw_experimental_setup: dict
) -> MeasurementSet:
    """Parse a measurements dataframe.

    :param measurement_table: result of running pd.read_csv on suitable file
    :param raw_experimental_setup: result of running toml.load on suitable file
    """
    experiments = [
        Experiment(**e) for e in raw_experimental_setup["experiment"]
    ]
    y = {
        mt: raw_measurement_table.loc[
            lambda df: df["measurement_type"] == mt.value  # noqa: B023
        ]
        for mt in MeasurementType
    }
    enzyme_knockouts = (
        [
            EnzymeKnockout(
                experiment_id=eko["experiment_id"], enzyme_id=eko["enzyme_id"]
            )
            for eko in raw_experimental_setup["enzyme_knockout"]
        ]
        if "enzyme_knockout" in raw_experimental_setup.keys()
        else None
    )
    pme_knockouts = (
        [
            PhosphorylationModifyingEnzymeKnockout(
                experiment_id=pko["experiment_id"], pme_id=pko["pme_id"]
            )
            for pko in raw_experimental_setup["phos_knockout"]
        ]
        if "phos_knockout" in raw_experimental_setup.keys()
        else None
    )
    return MeasurementSet(
        yconc=y[MeasurementType.MIC],
        yflux=y[MeasurementType.FLUX],
        yenz=y[MeasurementType.ENZYME],
        enzyme_knockouts=enzyme_knockouts,
        pme_knockouts=pme_knockouts,
        experiments=experiments,
    )
