from enum import Enum
from hashlib import md5
from typing import Any

from attrs import define, field

from classeq.core.domain.dtos.clade import ClasseqClade

from .._do_clade_adherence_test_for_single_sequence._dtos import AdherenceResult


class CladeAdherenceResultStatus(str, Enum):
    """The status of the clade adherence test."""

    INCONCLUSIVE = "inconclusive"
    NEXT_ITERATION = "next-iteration"
    CONCLUSIVE_INGROUP = "conclusive-ingroup"
    CONCLUSIVE_OUTGROUP = "conclusive-outgroup"
    MAX_RESOLUTION_REACHED = "max-resolution-reached"

    @classmethod
    def _missing_(cls, _: object) -> Any:
        return CladeAdherenceResultStatus.INCONCLUSIVE


@define
class CladeAdherenceResult:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    clade: ClasseqClade = field()
    ingroup_adherence_test: AdherenceResult = field()
    sister_adherence_test: AdherenceResult = field()

    # ? ------------------------------------------------------------------------
    # ? Life cycle hook methods
    # ? ------------------------------------------------------------------------

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, CladeAdherenceResult):
            return NotImplemented
        return self.__hash__() == other.__hash__()

    def __ne__(self, other: object) -> bool:
        # Not strictly necessary, but to avoid having both x==y and x!=y True at
        # the same time.
        if not isinstance(other, CladeAdherenceResult):
            return NotImplemented
        return not (self.__hash__() == other.__hash__())

    def __hash__(self) -> int:
        return int(
            md5(
                "-".join(
                    [
                        self.clade.__hash__().__str__(),
                        str(self.ingroup_adherence_test),
                        str(self.sister_adherence_test),
                    ]
                ).encode("utf-8")
            ).hexdigest(),
            16,
        )

    def __str__(self) -> str:
        return (
            "CladeAdherenceResult( "
            + f"clade={self.clade}, "
            + f"ingroup={self.ingroup_adherence_test.__str__()}, "
            + f"sister={self.sister_adherence_test.__str__()} )"
        )

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def to_dict(self, omit_children: bool = False) -> dict[str, Any]:
        return {
            "clade": self.clade.to_dict(omit_children=omit_children),
            "ingroup_adherence_test": self.ingroup_adherence_test.to_dict(),
            "sister_adherence_test": self.sister_adherence_test.to_dict(),
        }
