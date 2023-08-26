from __future__ import annotations

from enum import Enum

from attrs import define, field


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
    adherence_p_value: float = field(default=-999)

    def __str__(self) -> str:
        return (
            "AdherenceResult( "
            + f"match_kmers={self.match_kmers}, "
            + f"q={self.query_kmers_size}, "
            + f"s={self.subject_kmers_size}, "
            + f"status={self.status.name}, "
            + "adherence_p_value=%0.6f )" % self.adherence_p_value
        )
