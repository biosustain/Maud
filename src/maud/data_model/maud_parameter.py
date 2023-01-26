"""Definitions of Stan variables and StanVariableSet."""

from typing import List, Optional, Union

from pydantic import Field, root_validator, validator
from pydantic.dataclasses import dataclass

from maud.data_model.hardcoding import ID_SEPARATOR
from maud.data_model.id_component import IdComponent
from maud.data_model.maud_init import Init, Init1d, Init2d, InitAtomInput
from maud.data_model.prior import Prior
from maud.data_model.prior_input import IndPriorAtomInput, PriorMVNInput
from maud.getting_priors import (
    get_ind_prior_1d,
    get_ind_prior_2d,
    get_mvn_prior,
)


@dataclass
class MaudParameter:
    """A parameter in Maud's statistical model."""

    name: str
    shape_names: List[str]
    ids: List[List[str]]
    id_components: List[List[IdComponent]]
    split_ids: Optional[List[List[List[str]]]]
    non_negative: bool
    default_scale: float
    default_loc: float
    init_input: Optional[List[InitAtomInput]]
    prior_input: Optional[Union[PriorMVNInput, List[IndPriorAtomInput]]]
    prior_in_test_model: bool
    prior_in_train_model: bool
    prior: Prior = Field(init=False, exclude=True)
    inits: Init = Field(init=False, exclude=True)

    @validator("id_components")
    def id_components_have_same_length(cls, v):
        """Make sure that id_components contains lists with the same length."""
        first_length = len(v[0])
        assert all([len(x) == first_length for x in v])
        return v

    @root_validator(pre=False, skip_on_failure=True)
    def split_ids_exist_if_needed(cls, values):
        """Check split ids exist when there are non-trivial id components."""
        if any(len(idc) > 1 for idc in values["id_components"]):
            assert values["split_ids"] is not None
        return values


@dataclass(init=False)
class Km(MaudParameter):
    """Parameter representing a model's Michaelis constants."""

    def __init__(self, ids, split_ids, prior_input, init_input):
        self.name = "km"
        self.ids = ids
        self.split_ids = split_ids
        self.prior_input = prior_input
        self.init_input = init_input
        self.shape_names = ["N_km"]
        self.id_components = [
            [
                IdComponent.ENZYME,
                IdComponent.METABOLITE,
                IdComponent.COMPARTMENT,
            ]
        ]
        self.non_negative = True
        self.default_loc = -0.69  # roughly 0.5
        self.default_scale = 1
        self.prior_in_test_model = False
        self.prior_in_train_model = True
        self.mic_ids = [
            ID_SEPARATOR.join([met_id, compartment_id])
            for _, met_id, compartment_id in zip(*self.split_ids[0])
        ]
        self.prior = get_ind_prior_1d(
            self.prior_input,
            self.ids[0],
            self.id_components,
            self.non_negative,
            self.default_loc,
            self.default_scale,
        )
        self.inits = Init1d(
            self.ids,
            self.id_components,
            self.prior,
            self.init_input,
            self.non_negative,
        )

    @validator("split_ids")
    def split_ids_must_have_right_shape(cls, v):
        """Check that there are the right number of split ids."""
        assert v is not None, "Km split_ids are None"
        assert len(v) == 1, "Km split ids have wrong length"
        assert len(v[0]) == 4
        return v


@dataclass(init=False)
class Kcat(MaudParameter):
    """Stan variable representing a model's turnover numbers."""

    def __init__(self, ids, split_ids, prior_input, init_input):
        self.name = "kcat"
        self.ids = ids
        self.shape_names = ["N_enzyme_reaction"]
        self.id_components = [[IdComponent.ENZYME, IdComponent.REACTION]]
        self.split_ids = split_ids
        self.prior_input = prior_input
        self.init_input = init_input
        self.non_negative = True
        self.default_loc = -0.69  # roughly 0.5
        self.default_scale = 1
        self.prior_in_test_model = False
        self.prior_in_train_model = True
        self.prior = get_ind_prior_1d(
            self.prior_input,
            self.ids[0],
            self.id_components,
            self.non_negative,
            self.default_loc,
            self.default_scale,
        )
        self.inits = Init1d(
            self.ids,
            self.id_components,
            self.prior,
            self.init_input,
            self.non_negative,
        )


@dataclass(init=False)
class Ki(MaudParameter):
    """Stan variable representing a model's inhibition constants."""

    def __init__(self, ids, split_ids, prior_input, init_input):
        self.name = "ki"
        self.ids = ids
        self.split_ids = split_ids
        self.prior_input = prior_input
        self.init_input = init_input
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
        self.default_loc = -0.69  # roughly 0.5
        self.default_scale = 1
        self.prior_in_test_model = False
        self.prior_in_train_model = True
        self.mic_ids = [
            ID_SEPARATOR.join([met_id, compartment_id])
            for _, _, met_id, compartment_id in zip(*self.split_ids[0])
        ]
        self.er_ids = [
            enzyme_id
            for enzyme_id, reaction_id, _, _ in zip(*self.split_ids[0])
        ]
        self.prior = get_ind_prior_1d(
            self.prior_input,
            self.ids[0],
            self.id_components,
            self.non_negative,
            self.default_loc,
            self.default_scale,
        )
        self.inits = Init1d(
            self.ids,
            self.id_components,
            self.prior,
            self.init_input,
            self.non_negative,
        )


@dataclass(init=False)
class Dgf(MaudParameter):
    """Stan variable representing a model's standard formation energies."""

    def __init__(self, ids, split_ids, prior_input, init_input):
        self.name = "dgf"
        self.ids = ids
        self.split_ids = split_ids
        self.prior_input = prior_input
        self.init_input = init_input
        self.shape_names = ["N_metabolite"]
        self.id_components = [[IdComponent.METABOLITE]]
        self.non_negative = False
        self.default_loc = 0
        self.default_scale = 10
        self.prior_in_test_model = False
        self.prior_in_train_model = True
        self.prior = get_mvn_prior(
            self.prior_input,
            self.ids[0],
            self.id_components,
            self.non_negative,
            self.default_loc,
            self.default_scale,
        )
        self.inits = Init1d(
            self.ids,
            self.id_components,
            self.prior,
            self.init_input,
            self.non_negative,
        )


@dataclass(init=False)
class DissociationConstant(MaudParameter):
    """Stan variable representing a model's dissociation constants."""

    def __init__(self, ids, split_ids, prior_input, init_input):
        self.name = "dissociation_constant"
        self.ids = ids
        self.split_ids = split_ids
        self.prior_input = prior_input
        self.init_input = init_input
        self.shape_names = ["N_aa"]
        self.id_components = [
            [
                IdComponent.ENZYME,
                IdComponent.METABOLITE,
                IdComponent.COMPARTMENT,
                IdComponent.MODIFICATION_TYPE,
            ]
        ]
        self.non_negative = True
        self.mic_ids = [
            ID_SEPARATOR.join([met_id, compartment_id])
            for _, met_id, compartment_id, mt in zip(*self.split_ids[0])
        ]
        self.enzyme_ids = self.split_ids[0][0]
        self.default_loc = -0.69  # roughly 0.5
        self.default_scale = 1
        self.prior_in_test_model = False
        self.prior_in_train_model = True
        self.prior = get_ind_prior_1d(
            self.prior_input,
            self.ids[0],
            self.id_components,
            self.non_negative,
            self.default_loc,
            self.default_scale,
        )
        self.inits = Init1d(
            self.ids,
            self.id_components,
            self.prior,
            self.init_input,
            self.non_negative,
        )


@dataclass(init=False)
class TransferConstant(MaudParameter):
    """Stan variable representing a model's transfer constants."""

    def __init__(self, ids, split_ids, prior_input, init_input):
        self.name = "transfer_constant"
        self.ids = ids
        self.shape_names = ["N_ae"]
        self.id_components = [[IdComponent.ENZYME]]
        self.split_ids = split_ids
        self.prior_input = prior_input
        self.init_input = init_input
        self.non_negative = True
        self.default_loc = -0.69  # roughly 0.5
        self.default_scale = 1
        self.prior_in_test_model = False
        self.prior_in_train_model = True
        self.prior = get_ind_prior_1d(
            self.prior_input,
            self.ids[0],
            self.id_components,
            self.non_negative,
            self.default_loc,
            self.default_scale,
        )
        self.inits = Init1d(
            self.ids,
            self.id_components,
            self.prior,
            self.init_input,
            self.non_negative,
        )


@dataclass(init=False)
class KcatPme(MaudParameter):
    """Stan variable representing Kcats of phosphorylation modifying enzymes."""

    def __init__(self, ids, split_ids, prior_input, init_input):
        self.name = "kcat_pme"
        self.ids = ids
        self.shape_names = ["N_pme"]
        self.split_ids = split_ids
        self.prior_input = prior_input
        self.init_input = init_input
        self.id_components = [[IdComponent.PHOSPHORYLATION_MODIFYING_ENZYME]]
        self.non_negative = True
        self.default_loc = -0.69  # roughly 0.5
        self.default_scale = 1
        self.prior_in_test_model = False
        self.prior_in_train_model = True
        self.prior = get_ind_prior_1d(
            self.prior_input,
            self.ids[0],
            self.id_components,
            self.non_negative,
            self.default_loc,
            self.default_scale,
        )
        self.inits = Init1d(
            self.ids,
            self.id_components,
            self.prior,
            self.init_input,
            self.non_negative,
        )


@dataclass(init=False)
class Drain(MaudParameter):
    """Stan variable type for drain parameters."""

    def __init__(self, ids, split_ids, prior_input, init_input):
        self.ids = ids
        self.id_components = [[IdComponent.EXPERIMENT], [IdComponent.REACTION]]
        self.split_ids = split_ids
        self.prior_input = prior_input
        self.init_input = init_input
        self.non_negative = False
        self.default_loc = 0
        self.default_scale = 1
        self.prior = get_ind_prior_2d(
            self.prior_input,
            self.ids,
            self.id_components,
            self.non_negative,
            self.default_loc,
            self.default_scale,
        )
        self.inits = Init2d(
            self.ids,
            self.id_components,
            self.prior,
            self.init_input,
            self.non_negative,
        )


@dataclass(init=False)
class DrainTrain(Drain):
    """Stan variable for drain parameters of training experiments."""

    def __init__(self, ids, split_ids, prior_input, init_input):
        self.name = "drain_train"
        self.shape_names = ["N_experiment_train", "N_drain"]
        self.prior_in_test_model = False
        self.prior_in_train_model = True
        super().__init__(ids, split_ids, prior_input, init_input)


@dataclass(init=False)
class DrainTest(Drain):
    """Stan variable for drain parameters of test experiments."""

    def __init__(self, ids, split_ids, prior_input, init_input):
        self.name = "drain_test"
        self.shape_names = ["N_experiment_test", "N_drain"]
        self.prior_in_test_model = True
        self.prior_in_train_model = False
        super().__init__(ids, split_ids, prior_input, init_input)


@dataclass(init=False)
class ConcEnzyme(MaudParameter):
    """Parent class for enzyme concentration parameters."""

    def __init__(self, ids, split_ids, prior_input, init_input, measurements):
        self.ids = ids
        self.id_components = [[IdComponent.EXPERIMENT], [IdComponent.ENZYME]]
        self.split_ids = split_ids
        self.prior_input = prior_input
        self.init_input = init_input
        self.non_negative = True
        self.default_loc = -2.3
        self.default_scale = 2.0
        self.default_loc = 0.5
        self.default_scale = 1
        self.prior = get_ind_prior_2d(
            self.prior_input,
            self.ids,
            self.id_components,
            self.non_negative,
            self.default_loc,
            self.default_scale,
        )
        self.inits = Init2d(
            self.ids,
            self.id_components,
            self.prior,
            self.init_input,
            self.non_negative,
            measurements,
        )


@dataclass(init=False)
class ConcEnzymeTrain(ConcEnzyme):
    """Enzyme concentration parameters in training experiments."""

    def __init__(self, ids, split_ids, prior_input, init_input, measurements):
        self.name = "conc_enzyme_train"
        self.shape_names = ["N_experiment_train", "N_enzyme"]
        self.prior_in_test_model = False
        self.prior_in_train_model = True
        super().__init__(ids, split_ids, prior_input, init_input, measurements)


@dataclass(init=False)
class ConcEnzymeTest(ConcEnzyme):
    """Enzyme concentration parameters in test experiments."""

    def __init__(self, ids, split_ids, prior_input, init_input, measurements):
        self.name = "conc_enzyme_test"
        self.shape_names = ["N_experiment_test", "N_enzyme"]
        self.prior_in_test_model = True
        self.prior_in_train_model = False
        super().__init__(ids, split_ids, prior_input, init_input, measurements)


@dataclass(init=False)
class ConcUnbalanced(MaudParameter):
    """Parent class for unbalanced mic concentration parameters."""

    def __init__(self, ids, split_ids, prior_input, init_input, measurements):
        self.ids = ids
        self.split_ids = split_ids
        self.prior_input = prior_input
        self.init_input = init_input
        self.id_components = [
            [IdComponent.EXPERIMENT],
            [IdComponent.METABOLITE, IdComponent.COMPARTMENT],
        ]
        self.non_negative = True
        self.default_loc = -2.3
        self.default_scale = 2.0
        self.prior = get_ind_prior_2d(
            self.prior_input,
            self.ids,
            self.id_components,
            self.non_negative,
            self.default_loc,
            self.default_scale,
        )
        self.inits = Init2d(
            self.ids,
            self.id_components,
            self.prior,
            self.init_input,
            self.non_negative,
            measurements,
        )


@dataclass(init=False)
class ConcUnbalancedTrain(ConcUnbalanced):
    """Unbalanced mic concentration parameters in training experiments."""

    def __init__(self, ids, split_ids, prior_input, init_input, measurements):
        self.name = "conc_unbalanced_train"
        self.shape_names = ["N_experiment_train", "N_unbalanced"]
        self.prior_in_test_model = False
        self.prior_in_train_model = True
        super().__init__(ids, split_ids, prior_input, init_input, measurements)


@dataclass(init=False)
class ConcUnbalancedTest(ConcUnbalanced):
    """Unbalanced mic concentration parameters in test experiments."""

    def __init__(self, ids, split_ids, prior_input, init_input, measurements):
        self.name = "conc_unbalanced_test"
        self.shape_names = ["N_experiment_test", "N_enzyme"]
        self.prior_in_test_model = True
        self.prior_in_train_model = False
        super().__init__(ids, split_ids, prior_input, init_input, measurements)


@dataclass(init=False)
class ConcPme(MaudParameter):
    """Parent class for pme concentration parameters."""

    def __init__(self, ids, split_ids, prior_input, init_input):
        self.ids = ids
        self.split_ids = split_ids
        self.prior_input = prior_input
        self.init_input = init_input
        self.id_components = [
            [IdComponent.EXPERIMENT],
            [IdComponent.PHOSPHORYLATION_MODIFYING_ENZYME],
        ]
        self.non_negative = True
        self.default_loc = 0.1
        self.default_scale = 2.0
        self.prior = get_ind_prior_2d(
            self.prior_input,
            self.ids,
            self.id_components,
            self.non_negative,
            self.default_loc,
            self.default_scale,
        )
        self.inits = Init2d(
            self.ids,
            self.id_components,
            self.prior,
            self.init_input,
            self.non_negative,
        )


@dataclass(init=False)
class ConcPmeTrain(ConcPme):
    """Pme concentration parameters in training experiments."""

    def __init__(self, ids, split_ids, prior_input, init_input):
        self.name = "conc_pme_train"
        self.shape_names = ["N_experiment_train", "N_pme"]
        self.prior_in_test_model = False
        self.prior_in_train_model = True
        super().__init__(ids, split_ids, prior_input, init_input)


@dataclass(init=False)
class ConcPmeTest(ConcPme):
    """Pme concentration parameters in test experiments."""

    def __init__(self, ids, split_ids, prior_input, init_input):
        self.name = "conc_pme_test"
        self.shape_names = ["N_experiment_test", "N_pme"]
        self.prior_in_test_model = True
        self.prior_in_train_model = False
        super().__init__(ids, split_ids, prior_input, init_input)


@dataclass(init=False)
class Psi(MaudParameter):
    """Stan variable representing per-experiment membrane potentials."""

    def __init__(self, ids, split_ids, prior_input, init_input):
        self.ids = ids
        self.split_ids = split_ids
        self.prior_input = prior_input
        self.init_input = init_input
        self.id_components = [[IdComponent.EXPERIMENT]]
        self.non_negative = False
        self.default_loc = 0
        self.default_scale = 2
        self.prior = get_ind_prior_1d(
            self.prior_input,
            self.ids[0],
            self.id_components,
            self.non_negative,
            self.default_loc,
            self.default_scale,
        )
        self.inits = Init1d(
            self.ids,
            self.id_components,
            self.prior,
            self.init_input,
            self.non_negative,
        )


@dataclass(init=False)
class PsiTrain(Psi):
    """Pme concentration parameters in training experiments."""

    def __init__(self, ids, split_ids, prior_input, init_input):
        self.name = "psi_train"
        self.shape_names = ["N_experiment_train"]
        self.prior_in_test_model = False
        self.prior_in_train_model = True
        super().__init__(ids, split_ids, prior_input, init_input)


@dataclass(init=False)
class PsiTest(Psi):
    """Pme concentration parameters in test experiments."""

    def __init__(self, ids, split_ids, prior_input, init_input):
        self.name = "psi_test"
        self.shape_names = ["N_experiment_test"]
        self.prior_in_test_model = True
        self.prior_in_train_model = False
        super().__init__(ids, split_ids, prior_input, init_input)


@dataclass
class ParameterSet:
    """The parameters of a Maud input."""

    dgf: Dgf
    km: Km
    ki: Ki
    kcat: Kcat
    dissociation_constant: DissociationConstant
    transfer_constant: TransferConstant
    kcat_pme: KcatPme
    drain_train: DrainTrain
    drain_test: DrainTest
    conc_enzyme_train: ConcEnzymeTrain
    conc_enzyme_test: ConcEnzymeTest
    conc_unbalanced_train: ConcUnbalancedTrain
    conc_unbalanced_test: ConcUnbalancedTest
    conc_pme_train: ConcPmeTrain
    conc_pme_test: ConcPmeTest
    psi_train: PsiTrain
    psi_test: PsiTest
