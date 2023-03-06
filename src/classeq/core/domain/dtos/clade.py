from enum import Enum
from hashlib import md5
from typing import List, Self, Tuple
from uuid import UUID, uuid4

from attr import field, define

from classeq.core.domain.dtos.kmer_inverse_index import KmerIndex


class NodeType(Enum):
    ROOT = "root"
    OUTGROUP = "outgroup"
    INTERNAL = "internal"
    TERMINAL = "terminal"


@define(kw_only=True)
class CladeWrapper:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    name: str = field()
    id: UUID = field()
    type: NodeType = field()
    support: float | None = field(default=None)
    branch_length: float | None = field(default=None)
    parent: UUID | None = field(default=None)
    children: Tuple[Self, ...] | None = field(default=None)

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

    def __repr__(self) -> str:
        str_repr: List[str] = []

        backbone_color = "\033[90m"
        escape = "\033[0m"

        PIPE = f"{backbone_color}│{escape}"
        ELBOW = f"{backbone_color}└──{escape}"
        TEE = f"{backbone_color}├──{escape}"
        PIPE_PREFIX = f"{backbone_color}│   {escape}"
        SPACE_PREFIX = f"{backbone_color}    {escape}"

        def __print_nodes(clade: CladeWrapper, prefix: str) -> None:
            if clade.children is None:
                return

            entries = sorted(
                clade.children,
                key=lambda x: (
                    len(x.children) if x.children is not None else 1,
                    (
                        "".join([str(i.name) for i in x.children])
                        if x.children is not None
                        else "a"
                    ),
                    -x.support if x.support is not None else 0,
                    x.type.name,
                    x.name,
                ),
            )

            entries_length = len(entries)

            for index, child in enumerate(entries):
                if child.is_terminal() or child.is_outgroup():
                    base_color = "\033[92m"

                    if child.is_outgroup():
                        base_color = "\033[91m"

                    colored_name = f"{base_color}{child.name}{escape}"

                    if entries_length - 1 == index:
                        str_repr.append(prefix + ELBOW + f" {colored_name}")
                        str_repr.append(prefix + SPACE_PREFIX)
                        continue

                    str_repr.append(prefix + TEE + f" {colored_name}")
                    continue

                if child.is_internal():
                    indent_sep = PIPE_PREFIX
                    node_sep = TEE
                    int_support = (
                        None if child.support is None else round(child.support)
                    )

                    colored_support = f"\033[37m{int_support}{escape}"
                    colored_id = f"\033[94m{child.id}{escape}"

                    if entries_length - 1 == index:
                        indent_sep = SPACE_PREFIX
                        node_sep = ELBOW

                    str_repr.append(
                        prefix + node_sep + f"{colored_support} {colored_id}"
                    )

                    str_repr.append(prefix + indent_sep + PIPE)

                    __print_nodes(clade=child, prefix=prefix + indent_sep)

        str_repr.append(PIPE)
        __print_nodes(clade=self, prefix="")

        return "\n".join(str_repr)

    # ? ------------------------------------------------------------------------
    # ? Validations
    # ? ------------------------------------------------------------------------

    @id.default
    def _id_default(self) -> UUID:
        return uuid4()

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def is_root(self) -> bool:
        return self.type == NodeType.ROOT

    def is_outgroup(self) -> bool:
        return self.type == NodeType.OUTGROUP

    def is_internal(self) -> bool:
        return self.type == NodeType.INTERNAL

    def is_terminal(self) -> bool:
        return self.type == NodeType.TERMINAL
