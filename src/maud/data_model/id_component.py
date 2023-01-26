"""Provides Enum IdComponent defining things that can be ids in Maud."""

from enum import Enum


class IdComponent(str, Enum):
    """Possible kinds of id."""

    METABOLITE = "metabolite"
    COMPARTMENT = "compartment"
    ENZYME = "enzyme"
    REACTION = "reaction"
    EXPERIMENT = "experiment"
    PHOSPHORYLATION_MODIFYING_ENZYME = "phosphorylation_modifying_enzyme"
    MODIFICATION_TYPE = "modification_type"
