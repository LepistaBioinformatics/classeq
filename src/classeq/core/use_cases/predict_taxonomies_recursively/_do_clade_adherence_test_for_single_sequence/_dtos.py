from __future__ import annotations

from enum import Enum
from typing import Any

from attrs import define, field


class AdherenceStatus(Enum):
    SUCCESS = "success"
    NOT_ENOUGH_PRIORS = "not_enough_priors"
    UNDEFINED = "undefined"


@define(kw_only=True)
class AdherenceResult:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    match_kmers: int = field()
    query_kmers_size: int = field()
    subject_kmers_size: int = field()
    status: AdherenceStatus = field(default=AdherenceStatus.UNDEFINED)

    # ? ------------------------------------------------------------------------
    # ? Life cycle hook methods
    # ? ------------------------------------------------------------------------

    def __str__(self) -> str:
        return (
            "AdherenceResult( "
            + f"match_kmers={self.match_kmers}, "
            + f"q={self.query_kmers_size}, "
            + f"s={self.subject_kmers_size}, "
            + f"status={self.status.name} )"
        )

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return {
            "match_kmers": self.match_kmers,
            "query_kmers_size": self.query_kmers_size,
            "subject_kmers_size": self.subject_kmers_size,
            "status": self.status.name,
        }
