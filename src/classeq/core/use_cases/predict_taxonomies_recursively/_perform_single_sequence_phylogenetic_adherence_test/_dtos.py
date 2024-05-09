from enum import Enum
from typing import Any

from attrs import field, define

from classeq.settings import PREDICTION_DICT_KEY, STATUS_DICT_KEY

from .._perform_adherence_test_of_child_clades import (
    CladeAdherenceResult,
    CladeAdherenceResultStatus,
)


class AdherenceTestResultGroup(Enum):
    OUTGROUP = "outgroup"
    INGROUP = "ingroup"


@define(kw_only=True)
class PredictionStep:
    depth: int = field()
    result: CladeAdherenceResult = field()

    def to_dict(self, omit_children: bool = False) -> dict[str, Any]:
        return {
            "depth": self.depth,
            "result": self.result.to_dict(omit_children=omit_children),
        }


@define(kw_only=True)
class PredictionResult:
    status: CladeAdherenceResultStatus = field()
    path: list[PredictionStep] = field()

    def to_dict(self, omit_children: bool = False) -> dict[str, Any]:
        return {
            STATUS_DICT_KEY: self.status.value,
            PREDICTION_DICT_KEY: [
                step.to_dict(omit_children=omit_children) for step in self.path
            ],
        }
