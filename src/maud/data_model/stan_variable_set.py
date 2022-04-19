from dataclasses import field
from enum import Enum
from typing import List, Optional

from pydantic import validator
from pydantic.dataclasses import dataclass

from maud.data_model.hardcoding import ID_SEPARATOR


class IdComponent(Enum):
    METABOLITE = "metabolite"
    COMPARTMENT = "compartment"
    ENZYME = "enzyme"
    REACTION = "reaction"
    EXPERIMENT = "experiment"


class StanVariableType(Enum):
    PARAMETER = "parameter"
    GENERATED_QUANTITY = "generated_quantity"


@dataclass
class StanVariable:
    name: str
    var_type: StanVariableType = field(init=False)
    shape_names: List[str]
    ids: List[List[str]]
    id_components: List[List[IdComponent]]
    split_ids: Optional[List[List[List[str]]]]
    non_negative: bool
    default_scale: float
    default_loc: float

    def __post_init__(self):
        self.var_type = StanVariableType.PARAMETER

    @validator("id_components")
    def id_components_have_same_length(cls, v):
        first_length = len(v[0])
        assert all([len(x) == first_length for x in v])


@dataclass
class Km(StanVariable):
    def __init__(self, ids, split_ids):
        self.name = "km"
        self.var_type = StanVariableType.PARAMETER
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
        self.mic_ids: List[str] = field(init=False)
        self.er_ids: List[str] = field(init=False)
        self.default_loc = 0.5
        self.default_scale = 1

    def __post_init__(self):
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
        assert v is not None, "Km split_ids are None"
        assert len(v) == 1, "Km split ids have wrong length"
        assert len(v[0]) == 4
        return v


@dataclass(init=False)
class Kcat(StanVariable):
    def __init__(self, ids, split_ids):
        self.name = "kcat"
        self.var_type = StanVariableType.PARAMETER
        self.ids = ids
        self.shape_names = ["N_enzyme_reaction"]
        self.id_components = [[IdComponent.ENZYME, IdComponent.REACTION]]
        self.split_ids = split_ids
        self.non_negative = True
        self.default_loc = 0.5
        self.default_scale = 1


@dataclass(init=False)
class Ki(StanVariable):
    def __init__(self, ids, split_ids):
        self.name = "ki"
        self.var_type = StanVariableType.PARAMETER
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
        self.mic_ids: List[str] = field(init=False)
        self.er_ids: List[str] = field(init=False)
        self.default_loc = 0.5
        self.default_scale = 1

    def __post_init__(self):
        self.mic_ids = [
            ID_SEPARATOR.join([met_id, compartment_id])
            for _, _, met_id, compartment_id in self.split_ids[0]
        ]
        self.er_ids = [
            enzyme_id for enzyme_id, reaction_id, _, _ in self.split_ids[0]
        ]


@dataclass(init=False)
class Dgf(StanVariable):
    def __init__(self, ids):
        self.name = "dgf"
        self.var_type = StanVariableType.PARAMETER
        self.ids = ids
        self.split_ids = None
        self.shape_names = ["N_metabolite"]
        self.id_components = [[IdComponent.METABOLITE]]
        self.non_negative = False
        self.default_loc = 0
        self.default_scale = 10


@dataclass(init=False)
class DissociationConstant(StanVariable):
    def __init__(self, ids, split_ids):
        self.name = "dissociation_constant"
        self.var_type = StanVariableType.PARAMETER
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
        self.mic_ids: List[str] = field(init=False)
        self.enzyme_ids: List[str] = field(init=False)
        self.default_loc = 0.5
        self.default_scale = 1

    def __post_init__(self):
        self.mic_ids = [
            ID_SEPARATOR.join([met_id, compartment_id])
            for _, met_id, compartment_id in self.split_ids[0]
        ]
        self.enzyme_ids = [enzyme_id for enzyme_id, _, _ in self.split_ids[0]]


@dataclass(init=False)
class TransferConstant(StanVariable):
    def __init__(self, ids):
        self.name = "transfer_constant"
        self.var_type = StanVariableType.PARAMETER
        self.ids = ids
        self.shape_names = ["N_ae"]
        self.id_components = [[IdComponent.ENZYME]]
        self.split_ids = None
        self.non_negative = True
        self.default_loc = 0.5
        self.default_scale = 1


@dataclass(init=False)
class KcatPhos(StanVariable):
    def __init__(self, ids):
        self.name = "kcat_phos"
        self.var_type = StanVariableType.PARAMETER
        self.ids = ids
        self.shape_names = ["N_phosphorylation_enzymes"]
        self.split_ids = None
        self.id_components = [[IdComponent.ENZYME]]
        self.non_negative = True
        self.default_loc = 0.5
        self.default_scale = 1


@dataclass(init=False)
class Drain(StanVariable):
    def __init__(self, ids):
        self.name = "drain"
        self.var_type = StanVariableType.PARAMETER
        self.ids = ids
        self.shape_names = ["N_experiment", "N_drain"]
        self.id_components = [[IdComponent.EXPERIMENT], [IdComponent.REACTION]]
        self.split_ids = None
        self.non_negative = False
        self.default_loc = 0
        self.default_scale = 1


@dataclass(init=False)
class ConcEnzyme(StanVariable):
    def __init__(self, ids):
        self.name = "conc_enzyme"
        self.var_type = StanVariableType.PARAMETER
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
    def __init__(self, ids):
        self.name = "conc_unbalanced"
        self.var_type = StanVariableType.PARAMETER
        self.ids = ids
        self.shape_names = ["N_experiment", "N_unbalanced"]
        self.split_ids = None
        self.id_components = [
            [
                IdComponent.EXPERIMENT,
                IdComponent.METABOLITE,
                IdComponent.COMPARTMENT,
            ]
        ]
        self.non_negative = True
        self.default_loc = 0.1
        self.default_scale = 2.0


@dataclass(init=False)
class ConcPhos(StanVariable):
    def __init__(self, ids):
        self.name = "conc_phos"
        self.var_type = StanVariableType.PARAMETER
        self.ids = ids
        self.shape_names = ["N_experiment", "N_phosphorylation"]
        self.split_ids = None
        self.id_components = [[IdComponent.EXPERIMENT, IdComponent.ENZYME]]
        self.non_negative = True
        self.default_loc = 0.1
        self.default_scale = 2.0


@dataclass
class StanVariableSet:
    dgf: Dgf
    km: Km
    ki: Ki
    kcat: Kcat
    dissociation_constant: DissociationConstant
    transfer_constant: TransferConstant
    kcat_phos: KcatPhos
    drain: Drain
    conc_enzyme: ConcEnzyme
    conc_unbalanced: ConcUnbalanced
    conc_phos: ConcPhos
