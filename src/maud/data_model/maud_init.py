from typing import Dict, Sequence, Union

InitValue = Union[
    float,
    Sequence[float],
    Sequence[Sequence[float]],
    Sequence[Sequence[Sequence[float]]]
]

InitDict = Dict[str, InitValue]
