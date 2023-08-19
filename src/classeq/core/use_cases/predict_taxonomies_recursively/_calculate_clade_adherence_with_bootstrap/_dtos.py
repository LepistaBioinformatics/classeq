from __future__ import annotations

from enum import Enum
from typing import Any

from attrs import define, field


class AdherenceTestStrategy(Enum):
    KMERS_INTERSECTION = "kmers-intersection"
    JOINT_PROBABILITY = "joint-probability"

    @classmethod
    def _missing_(cls, _: Any = None) -> AdherenceTestStrategy:
        return cls.JOINT_PROBABILITY


class AdherenceStatus(Enum):
    SUCCESS = "success"
    NOT_ENOUGH_PRIORS = "not_enough_priors"
    UNDEFINED = "undefined"


@define(kw_only=True)
class AdherenceResult:
    match_kmers: int = field()
    query_kmers_size: int = field()
    subject_kmers_size: int = field()
    status: AdherenceStatus = field(default=AdherenceStatus.UNDEFINED)
    joint_probability: float = field(default=-999)
    adherence_p_value: float = field(default=-999)
    strategy: AdherenceTestStrategy = field(
        default=AdherenceTestStrategy.JOINT_PROBABILITY
    )

    def __str__(self) -> str:
        return (
            "AdherenceResult( "
            + f"match_kmers={self.match_kmers}, "
            + f"q={self.query_kmers_size}, "
            + f"s={self.subject_kmers_size}, "
            + f"status={self.status.name}, "
            + "joint_probability=%0.10f, " % self.joint_probability
            + "adherence_p_value=%0.6f )" % self.adherence_p_value
        )
