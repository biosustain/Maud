"""Provides dataclass StanCoordSet"""
from typing import List

from pydantic.dataclasses import dataclass


@dataclass
class StanCoordSet:
    """Object containing human-readable indexes for Maud's parameters.

    These are "coordinates" in the sense of xarray

    """

    metabolite: List[str]
    mic: List[str]
    balanced_mic: List[str]
    unbalanced_mic: List[str]
    km: List[str]
    reaction: List[str]
    experiment: List[str]
    enzyme: List[str]
    edge: List[str]
    allosteric_enzyme: List[str]
    drain: List[str]
    phos_enz: List[str]
    yconc_exp: List[str]
    yconc_mic: List[str]
    yflux_exp: List[str]
    yflux_rxn: List[str]
    yenz_exp: List[str]
    yenz_enz: List[str]
    ci_enzs: List[str]
    ci_mics: List[str]
    ai_enzs: List[str]
    ai_mics: List[str]
    aa_enzs: List[str]
    aa_mics: List[str]
    enz_ko_exps: List[str]
    enz_ko_enzs: List[str]
    phos_ko_exps: List[str]
    phos_ko_enzs: List[str]
