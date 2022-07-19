"""Definitions of Stan variables and StanVariableSet."""

from enum import Enum
from typing import List, Optional

from pydantic import Field, root_validator, validator
from pydantic.dataclasses import dataclass

from maud.data_model.hardcoding import ID_SEPARATOR


class IdComponent(str, Enum):
    """Possible kinds of id."""

    METABOLITE = "metabolite"
    COMPARTMENT = "compartment"
    ENZYME = "enzyme"
    REACTION = "reaction"
    EXPERIMENT = "experiment"
    PHOSPHORYLATION_MODIFYING_ENZYME = "phosphorylation_modifying_enzyme"


@dataclass
class StanVariable:
    """A thing that could be a random variable in Maud's statistical model.

    Probably appearing the parameters, transformed parameters or generated
    quantities block of the Stan program.

    """

    name: str
    shape_names: List[str]
    ids: List[List[str]]
    id_components: List[List[IdComponent]]
    split_ids: Optional[List[List[List[str]]]]
    non_negative: bool
    default_scale: float
    default_loc: float

    @validator("id_components")
    def id_components_have_same_length(cls, v):
        """Make sure that id_components contains lists with the same length."""
        first_length = len(v[0])
        assert all([len(x) == first_length for x in v])
        return v

    @root_validator
    def split_ids_exist_if_needed(cls, values):
        """Check split ids exist when there are non-trivial id components."""
        if any(len(idc) > 1 for idc in values["id_components"]):
            assert values["split_ids"] is not None
        return values


@dataclass
class Km(StanVariable):
    """Stan variable representing a model's Michaelis constants."""

    def __init__(self, ids, split_ids):
        self.name = "km"
        self.ids = ids
        self.split_ids = split_ids
        self.shape_names = ["N_km"]
        self.id_components = [
            [
                IdComponent.ENZYME,
                IdComponent.METABOLITE,
                IdComponent.COMPARTMENT,
            ]
        ]
        self.non_negative = True
        self.mic_ids: List[str] = Field(init=False, exclude=True)
        self.er_ids: List[str] = Field(init=False, exclude=True)
        self.default_loc = 0.5
        self.default_scale = 1

    def __post_init__(self):
        """Set mic_ids and er_ids."""
        self.mic_ids = [
            ID_SEPARATOR.join([met_id, compartment_id])
            for _, met_id, compartment_id in self.split_ids[0]
        ]
        self.er_ids = [
            ID_SEPARATOR.join([enzyme_id, reaction_id])
            for enzyme_id, reaction_id, _, _ in self.split_ids[0]
        ]

    @validator("split_ids")
    def split_ids_must_have_right_shape(cls, v):
        """Check that there are the right number of split ids."""
        assert v is not None, "Km split_ids are None"
        assert len(v) == 1, "Km split ids have wrong length"
        assert len(v[0]) == 4
        return v


@dataclass(init=False)
class Kcat(StanVariable):
    """Stan variable representing a model's turnover numbers."""

    def __init__(self, ids, split_ids):
        self.name = "kcat"
        self.ids = ids
        self.shape_names = ["N_enzyme_reaction"]
        self.id_components = [[IdComponent.ENZYME, IdComponent.REACTION]]
        self.split_ids = split_ids
        self.non_negative = True
        self.default_loc = 0.5
        self.default_scale = 1


@dataclass(init=False)
class Ki(StanVariable):
    """Stan variable representing a model's inhibition constants."""

    def __init__(self, ids, split_ids):
        self.name = "ki"
        self.ids = ids
        self.split_ids = split_ids
        self.shape_names = ["N_ci"]
        self.id_components = [
            [
                IdComponent.ENZYME,
                IdComponent.REACTION,
                IdComponent.METABOLITE,
                IdComponent.COMPARTMENT,
            ]
        ]
        self.non_negative = True
        self.mic_ids: List[str] = Field(init=False, exclude=True)
        self.er_ids: List[str] = Field(init=False, exclude=True)
        self.default_loc = 0.5
        self.default_scale = 1

    def __post_init__(self):
        """Set mic ids and er ids."""
        self.mic_ids = [
            ID_SEPARATOR.join([met_id, compartment_id])
            for _, _, met_id, compartment_id in self.split_ids[0]
        ]
        self.er_ids = [
            enzyme_id for enzyme_id, reaction_id, _, _ in self.split_ids[0]
        ]


@dataclass(init=False)
class Dgf(StanVariable):
    """Stan variable representing a model's standard formation energies."""

    def __init__(self, ids):
        self.name = "dgf"
        self.ids = ids
        self.split_ids = None
        self.shape_names = ["N_metabolite"]
        self.id_components = [[IdComponent.METABOLITE]]
        self.non_negative = False
        self.default_loc = 0
        self.default_scale = 10


@dataclass(init=False)
class DissociationConstant(StanVariable):
    """Stan variable representing a model's dissociation constants."""

    def __init__(self, ids, split_ids):
        self.name = "dissociation_constant"
        self.ids = ids
        self.split_ids = split_ids
        self.shape_names = ["N_aa"]
        self.id_components = [
            [
                IdComponent.ENZYME,
                IdComponent.METABOLITE,
                IdComponent.COMPARTMENT,
            ]
        ]
        self.non_negative = True
        self.mic_ids: List[str] = Field(init=False, exclude=True)
        self.enzyme_ids: List[str] = Field(init=False, exclude=True)
        self.default_loc = 0.5
        self.default_scale = 1

    def __post_init__(self):
        """Set mic ids and enzyme ids."""
        self.mic_ids = [
            ID_SEPARATOR.join([met_id, compartment_id])
            for _, met_id, compartment_id in self.split_ids[0]
        ]
        self.enzyme_ids = [enzyme_id for enzyme_id, _, _ in self.split_ids[0]]


@dataclass(init=False)
class TransferConstant(StanVariable):
    """Stan variable representing a model's transfer constants."""

    def __init__(self, ids):
        self.name = "transfer_constant"
        self.ids = ids
        self.shape_names = ["N_ae"]
        self.id_components = [[IdComponent.ENZYME]]
        self.split_ids = None
        self.non_negative = True
        self.default_loc = 0.5
        self.default_scale = 1


@dataclass(init=False)
class KcatPme(StanVariable):
    """Stan variable representing Kcats of phosphorylation modifying enzymes."""

    def __init__(self, ids):
        self.name = "kcat_pme"
        self.ids = ids
        self.shape_names = ["N_pme"]
        self.split_ids = None
        self.id_components = [[IdComponent.PHOSPHORYLATION_MODIFYING_ENZYME]]
        self.non_negative = True
        self.default_loc = 0.5
        self.default_scale = 1


@dataclass(init=False)
class Drain(StanVariable):
    """Stan variable representing a model's drain parameters."""

    def __init__(self, ids):
        self.name = "drain"
        self.ids = ids
        self.shape_names = ["N_experiment", "N_drain"]
        self.id_components = [[IdComponent.EXPERIMENT], [IdComponent.REACTION]]
        self.split_ids = None
        self.non_negative = False
        self.default_loc = 0
        self.default_scale = 1


@dataclass(init=False)
class ConcEnzyme(StanVariable):
    """Stan variable representing a model's enzyme concentrations."""

    def __init__(self, ids):
        self.name = "conc_enzyme"
        self.ids = ids
        self.shape_names = ["N_experimeent", "N_enzyme"]
        self.id_components = [[IdComponent.EXPERIMENT], [IdComponent.ENZYME]]
        self.split_ids = None
        self.non_negative = True
        self.default_loc = 0.1
        self.default_scale = 2.0
        self.default_loc = 0.5
        self.default_scale = 1


@dataclass(init=False)
class ConcUnbalanced(StanVariable):
    """Stan variable representing a model's unbalanced mic concentrations."""

    def __init__(self, ids, split_ids):
        self.name = "conc_unbalanced"
        self.ids = ids
        self.shape_names = ["N_experiment", "N_unbalanced"]
        self.split_ids = split_ids
        self.id_components = [
            [IdComponent.EXPERIMENT],
            [IdComponent.METABOLITE, IdComponent.COMPARTMENT],
        ]
        self.non_negative = True
        self.default_loc = 0.1
        self.default_scale = 2.0


@dataclass(init=False)
class ConcPme(StanVariable):
    """Stan variable representing a model's pme concentrations."""

    def __init__(self, ids):
        self.name = "conc_pme"
        self.ids = ids
        self.shape_names = ["N_experiment", "N_pme"]
        self.split_ids = None
        self.id_components = [
            [IdComponent.EXPERIMENT],
            [IdComponent.PHOSPHORYLATION_MODIFYING_ENZYME],
        ]
        self.non_negative = True
        self.default_loc = 0.1
        self.default_scale = 2.0


@dataclass(init=False)
class Psi(StanVariable):
    """Stan variable representing per-experiment membrane potentials."""

    def __init__(self, ids):
        self.name = "psi"
        self.ids = ids
        self.shape_names = ["N_experiment"]
        self.split_ids = None
        self.id_components = [[IdComponent.EXPERIMENT]]
        self.non_negative = False
        self.default_loc = 0
        self.default_scale = 2


@dataclass
class StanVariableSet:
    """A MaudInput's set of Stan variables."""

    dgf: Dgf
    km: Km
    ki: Ki
    kcat: Kcat
    dissociation_constant: DissociationConstant
    transfer_constant: TransferConstant
    kcat_pme: KcatPme
    drain: Drain
    conc_enzyme: ConcEnzyme
    conc_unbalanced: ConcUnbalanced
    conc_pme: ConcPme
    psi: Psi
