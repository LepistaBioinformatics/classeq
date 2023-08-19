from enum import Enum
from hashlib import md5

from attrs import define, field

from classeq.core.domain.dtos.clade import ClasseqClade

from .._calculate_clade_adherence_with_bootstrap._dtos import AdherenceResult


class CladeAdherenceResultStatus(Enum):
    """The status of the clade adherence test.

    Description:
        This class defines the status of the clade adherence test. The status
        can be one of the following:
            - `next-iteration`: The test is inconclusive and the test should
                continue to the next iteration.
            - `conclusive-ingroup`: The test is conclusive and the sequence is
                in the ingroup.
            - `conclusive-sister`: The test is conclusive and the sequence is
                in the sister group.
            - `inconclusive`: The test is inconclusive and the sequence is not
                in the ingroup or sister group.
            - `not-possible`: The test is not possible because the clade is a
                leaf node.

    Attributes:
        NEXT_ITERATION (str): The test is inconclusive and the test should
            continue to the next iteration.
        CONCLUSIVE_INGROUP (str): The test is conclusive and the sequence is
            in the ingroup.
        CONCLUSIVE_SISTER (str): The test is conclusive and the sequence is
            in the sister group.
        MAX_DEPTH_REACHED (str): The maximum depth of the tree was reached.
        NOT_POSSIBLE (str): The test is not possible because the clade is a
            leaf node.

    """

    INCONCLUSIVE = "inconclusive"
    NEXT_ITERATION = "next-iteration"
    CONCLUSIVE_INGROUP = "conclusive-ingroup"
    MAX_RESOLUTION_REACHED = "max-resolution-reached"


@define
class CladeAdherenceResult:
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