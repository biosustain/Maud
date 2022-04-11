"""Provides dataclass SimulationStudyOutput"""
from cmdstanpy import CmdStanMCMC
from pydantic.dataclasses import dataclass

from maud.data_model.maud_input import MaudInput


@dataclass
class SimulationStudyOutput:
    """Expected output of a simulation study.

    :param input_data_sim: dictionary used to create simulation
    :param input_data_sample: dictionary used to create samples
    :param true_values: dictionary mapping param names to true values
    :param sim: CmdStanMCMC that generated simulated measurements
    :param mi: Maud input used for sampling
    :param samples: CmdStanMCMC output of the simulation study
    """

    def __init__(
        self,
        input_data_sim: dict,
        input_data_sample: dict,
        true_values: dict,
        sim: CmdStanMCMC,
        mi: MaudInput,
        samples: CmdStanMCMC,
    ):
        self.input_data_sim = input_data_sim
        self.input_data_sample = input_data_sample
        self.true_values = true_values
        self.sim = sim
        self.mi = mi
        self.samples = samples
