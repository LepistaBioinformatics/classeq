from enum import Enum
from hashlib import md5
from uuid import UUID, uuid4

from attr import field, frozen

from classeq.core.domain.dtos.kmer_inverse_index import KmerIndex


class NodeType(Enum):
    ROOT = "root"
    OUTGROUP = "outgroup"
    INTERNAL = "internal"
    TERMINAL = "terminal"


@frozen(kw_only=True)
class CladeWrapper:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    name: str = field()
    id: UUID = field()
    type: NodeType = field()
    support: float | None = field(default=None)
    parent: UUID | None = field(default=None)

    # ? ------------------------------------------------------------------------
    # ? Life cycle hook methods
    # ? ------------------------------------------------------------------------

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, KmerIndex):
            return NotImplemented
        return self.__hash__() == other.__hash__()

    def __ne__(self, other: object) -> bool:
        # Not strictly necessary, but to avoid having both x==y and x!=y True at
        # the same time.
        return not (self.__hash__() == other.__hash__())

    def __hash__(self) -> int:
        return int(md5(self.id.__str__().encode("utf-8")).hexdigest(), 16)

    # ? ------------------------------------------------------------------------
    # ? Validations
    # ? ------------------------------------------------------------------------

    @id.default
    def _id_default(self) -> UUID:
        return uuid4()
