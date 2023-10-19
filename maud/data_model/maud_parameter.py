"""Provides model MaudParameter, and subclasses for all parameters Maud uses."""

from copy import deepcopy
from typing import List, Optional, Union

from pydantic import BaseModel, computed_field, field_validator, model_validator

from maud.data_model.experiment import Measurement
from maud.data_model.hardcoding import ID_SEPARATOR
from maud.data_model.id_component import IdComponent
from maud.data_model.maud_init import Init, Init1d, Init2d, InitAtomInput
from maud.data_model.parameter_input import (
    ParameterInputAtom,
    ParameterInputMVN,
)
from maud.data_model.prior import IndPrior1d, IndPrior2d, PriorMVN


class MaudParameter(BaseModel):
    """A parameter in Maud's statistical model."""

    name: str
    shape_names: List[str]
    id_components: List[List[IdComponent]]
    non_negative: bool
    default_scale: float
    default_loc: float
    prior_in_test_model: bool
    prior_in_train_model: bool
    user_input: Optional[Union[List[ParameterInputAtom], ParameterInputMVN]]
    init_input: Optional[List[InitAtomInput]]
    ids: List[List[str]]
    split_ids: List[List[List[str]]]
    measurements: Optional[List[Measurement]] = None

    @computed_field
    def fixed_ids(self) -> Optional[List[List[str]]]:
        """Set the fixed_ids field."""
        if self.name != "dgf":
            return None
        elif self.user_input is None:
            return None
        elif isinstance(self.user_input, List):
            out = [[]]
            for pia in self.user_input:
                if pia.fixed_value is not None:
                    out[0] += [
                        ID_SEPARATOR.join([getattr(pia, c) for c in idci])
                        for idci in self.id_components
                    ]
            return out
        elif isinstance(self.user_input, ParameterInputMVN):
            if self.user_input.fixed_values is None:
                return None
            else:
                return [list(self.user_input.fixed_values.keys())]
        else:
            raise ValueError(f"Something wrong with input {self.user_input}")

    @computed_field
    def fixed_values(self) -> Optional[List[List[float]]]:
        """Set the fixed_values field."""
        if self.name != "dgf":
            return None
        elif self.user_input is None:
            return None
        elif isinstance(self.user_input, List):
            out = [[]]
            for pia in self.user_input:
                if pia.fixed_value is not None:
                    out[0].append(pia.fixed_value)
            return out
        elif isinstance(self.user_input, ParameterInputMVN):
            if self.user_input.fixed_values is None:
                return None
            else:
                return [list(self.user_input.fixed_values.values())]
        else:
            raise ValueError(f"Something wrong with input {self.user_input}")

    @computed_field
    def prior(self) -> Union[IndPrior1d, IndPrior2d, PriorMVN]:
        """Return a prior, calculated from the user input."""
        if self.name == "dgf":
            initialiser = PriorMVN
        elif len(self.shape_names) == 1:
            initialiser = IndPrior1d
        else:
            initialiser = IndPrior2d
        return initialiser(
            user_input=self.user_input,
            ids=self.ids,
            id_components=self.id_components,
            non_negative=self.non_negative,
            default_loc=self.default_loc,
            default_scale=self.default_scale,
        )

    @computed_field
    def inits(self) -> Init:
        """Add the inits field."""
        initialiser = Init1d if len(self.shape_names) == 1 else Init2d
        return initialiser(
            ids=self.ids,
            id_components=self.id_components,
            prior=self.prior,
            param_init_input=self.init_input,
            non_negative=self.non_negative,
            measurements=self.measurements,
        )

    @field_validator("id_components")
    def id_components_have_same_length(cls, v):
        """Make sure that id_components contains lists with the same length."""
        first_length = len(v[0])
        assert all([len(x) == first_length for x in v])
        return v

    @model_validator(mode="after")
    def split_ids_exist_if_needed(self):
        """Check split ids exist when there are non-trivial id components."""
        if any(len(idc) > 1 for idc in self.id_components):
            assert self.split_ids is not None
        return self


class TrainTestParameter(MaudParameter):
    """Mark parameter to have different priors between train and test.

    The class must be filled with `shape_names` and `name` for train.
    The test version is created by calling the `test` method.
    """

    prior_in_test_model: bool = False
    prior_in_train_model: bool = True

    def test(self):
        """Generate the test counterpart."""
        test_self = deepcopy(self)
        test_self.shape_names = [
            shape_name.replace("train", "test")
            for shape_name in self.shape_names
        ]
        test_self.name = self.name.replace("train", "test")
        test_self.prior_in_test_model = True
        test_self.prior_in_train_model = False
        return test_self


class Km(MaudParameter):
    """Parameter representing a model's Michaelis constants."""

    name: str = "km"
    shape_names: List[str] = ["N_km"]
    id_components: List[List[IdComponent]] = [
        [
            IdComponent.ENZYME,
            IdComponent.METABOLITE,
            IdComponent.COMPARTMENT,
        ]
    ]
    non_negative: bool = True
    default_loc: float = -0.69  # roughly 0.5
    default_scale: float = 1
    prior_in_test_model: bool = False
    prior_in_train_model: bool = True

    @computed_field
    def mic_ids(self) -> List[str]:
        """Add the mic_ids field."""
        return [
            ID_SEPARATOR.join([met_id, compartment_id])
            for _, met_id, compartment_id in zip(*self.split_ids[0])
        ]

    @field_validator("split_ids")
    def split_ids_must_have_right_shape(cls, v):
        """Check that there are the right number of split ids."""
        assert v is not None, "Km split_ids are None"
        assert len(v) == 1, "Km split ids have wrong length"
        assert len(v[0]) == 3, f"Wrong number of split id components: {v[0]}"
        return v


class Kcat(MaudParameter):
    """Stan variable representing a model's turnover numbers."""

    name: str = "kcat"
    shape_names: List[str] = ["N_enzyme_reaction"]
    id_components: List[List[IdComponent]] = [
        [IdComponent.ENZYME, IdComponent.REACTION]
    ]
    non_negative: bool = True
    default_loc: float = -0.69  # roughly 0.5
    default_scale: float = 1
    prior_in_test_model: bool = False
    prior_in_train_model: bool = True


class Ki(MaudParameter):
    """Stan variable representing a model's inhibition constants."""

    name: str = "ki"
    shape_names: List[str] = ["N_ci"]
    id_components: List[List[IdComponent]] = [
        [
            IdComponent.ENZYME,
            IdComponent.REACTION,
            IdComponent.METABOLITE,
            IdComponent.COMPARTMENT,
        ]
    ]
    non_negative: bool = True
    default_loc: float = -0.69  # roughly 0.5
    default_scale: float = 1
    prior_in_test_model: bool = False
    prior_in_train_model: bool = True

    @computed_field
    def mic_ids(self) -> List[str]:
        """Get the mic ids."""
        return [
            ID_SEPARATOR.join([met_id, compartment_id])
            for _, _, met_id, compartment_id in zip(*self.split_ids[0])
        ]

    @computed_field
    def er_ids(self) -> List[str]:
        """Get the enzyme-reaction ids."""
        return [
            ID_SEPARATOR.join([enzyme_id, reaction_id])
            for enzyme_id, reaction_id, _, _ in zip(*self.split_ids[0])
        ]


class Dgf(MaudParameter):
    """Stan variable representing a model's standard formation energies."""

    name: str = "dgf"
    shape_names: List[str] = ["N_metabolite"]
    id_components: List[List[IdComponent]] = [[IdComponent.METABOLITE]]
    non_negative: bool = False
    default_loc: float = 0
    default_scale: float = 10
    prior_in_test_model: bool = False
    prior_in_train_model: bool = True


class DissociationConstant(MaudParameter):
    """Stan variable representing a model's dissociation constants."""

    name: str = "dissociation_constant"
    shape_names: List[str] = ["N_aa"]
    id_components: List[List[IdComponent]] = [
        [
            IdComponent.ENZYME,
            IdComponent.METABOLITE,
            IdComponent.COMPARTMENT,
            IdComponent.MODIFICATION_TYPE,
        ]
    ]
    non_negative: bool = True
    default_loc: float = -0.69  # roughly 0.5
    default_scale: float = 1
    prior_in_test_model: bool = False
    prior_in_train_model: bool = True

    @computed_field
    def mic_ids(self) -> List[str]:
        """Get the mic ids."""
        return [
            ID_SEPARATOR.join([met_id, compartment_id])
            for _, _, met_id, compartment_id in zip(*self.split_ids[0])
        ]

    @computed_field
    def enzyme_ids(self) -> List[str]:
        """Get the enzyme ids."""
        return self.split_ids[0][0]


class TransferConstant(MaudParameter):
    """Stan variable representing a model's transfer constants."""

    name: str = "transfer_constant"
    shape_names: List[str] = ["N_ae"]
    id_components: List[List[IdComponent]] = [[IdComponent.ENZYME]]
    non_negative: bool = True
    default_loc: float = -0.69  # roughly 0.5
    default_scale: float = 1
    prior_in_test_model: bool = False
    prior_in_train_model: bool = True


class KcatPme(MaudParameter):
    """Stan variable representing Kcats of phosphorylation modifying enzymes."""

    name: str = "kcat_pme"
    shape_names: List[str] = ["N_pme"]
    id_components: List[List[IdComponent]] = [
        [IdComponent.PHOSPHORYLATION_MODIFYING_ENZYME]
    ]
    non_negative: bool = True
    default_loc: float = -0.69  # roughly 0.5
    default_scale: float = 1
    prior_in_test_model: bool = False
    prior_in_train_model: bool = True


class Drain(TrainTestParameter):
    """Stan variable type for drain parameters."""

    name: str = "drain_train"
    shape_names: List[str] = ["N_experiment_train", "N_drain"]
    id_components: List[List[IdComponent]] = [
        [IdComponent.EXPERIMENT],
        [IdComponent.REACTION],
    ]
    non_negative: bool = False
    default_loc: float = 0
    default_scale: float = 1


class ConcEnzyme(TrainTestParameter):
    """Parent class for enzyme concentration parameters."""

    name: str = "conc_enzyme_train"
    shape_names: List[str] = ["N_experiment_train", "N_enzyme"]
    id_components: List[List[IdComponent]] = [
        [IdComponent.EXPERIMENT],
        [IdComponent.ENZYME],
    ]
    non_negative: bool = True
    default_loc: float = -2.3
    default_scale: float = 2.0
    default_loc: float = 0.5
    default_scale: float = 1


class ConcUnbalanced(TrainTestParameter):
    """Parent class for unbalanced mic concentration parameters."""

    name: str = "conc_unbalanced_train"
    shape_names: List[str] = ["N_experiment_train", "N_unbalanced"]
    id_components: List[List[IdComponent]] = [
        [IdComponent.EXPERIMENT],
        [IdComponent.METABOLITE, IdComponent.COMPARTMENT],
    ]
    non_negative: bool = True
    default_loc: float = -2.3
    default_scale: float = 2.0


class ConcPme(TrainTestParameter):
    """Parent class for pme concentration parameters."""

    name: str = "conc_pme_train"
    shape_names: List[str] = ["N_experiment_train", "N_pme"]
    id_components: List[List[IdComponent]] = [
        [IdComponent.EXPERIMENT],
        [IdComponent.PHOSPHORYLATION_MODIFYING_ENZYME],
    ]
    non_negative: bool = True
    default_loc: float = 0.1
    default_scale: float = 2.0


class Psi(TrainTestParameter):
    """Stan variable representing per-experiment membrane potentials."""

    name: str = "psi_train"
    shape_names: List[str] = ["N_experiment_train"]
    id_components: List[List[IdComponent]] = [[IdComponent.EXPERIMENT]]
    non_negative: bool = False
    default_loc: float = 0
    default_scale: float = 2
