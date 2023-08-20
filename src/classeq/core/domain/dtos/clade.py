from enum import Enum
from hashlib import md5
from typing import Any, Self
from uuid import UUID, uuid4

import clean_base.exceptions as c_exc
from attr import define, field
from clean_base.either import Either, left, right

from classeq.core.domain.dtos.kmer_inverse_index import KmerIndex
from classeq.settings import LOGGER


class NodeType(Enum):
    ROOT = "root"
    OUTGROUP = "outgroup"
    INTERNAL = "internal"
    TERMINAL = "terminal"


@define(kw_only=True)
class ClasseqClade:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    name: str = field()
    id: UUID = field()
    type: NodeType = field()
    support: float | None = field(default=None)
    branch_length: float | None = field(default=None)
    parent: UUID | None = field(default=None)
    children: tuple[Self, ...] | None = field(default=None)

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
        return self.get_pretty_clade()

    # ? ------------------------------------------------------------------------
    # ? Validations
    # ? ------------------------------------------------------------------------

    @id.default
    def _id_default(self) -> UUID:
        return uuid4()

    # ? ------------------------------------------------------------------------
    # ? Public class methods
    # ? ------------------------------------------------------------------------

    @classmethod
    def from_dict(
        cls,
        content: dict[str, Any],
    ) -> Either[c_exc.MappedErrors, Self]:
        for key in [
            "name",
            "id",
            "type",
            "support",
            "branch_length",
            "parent",
            "children",
        ]:
            if key not in content:
                return left(
                    c_exc.DadaTransferObjectError(
                        f"Invalid content detected on parse `{ClasseqClade}`."
                        f" {key}` key is empty.",
                        logger=LOGGER,
                    )
                )

        children: list[Self] = []
        if isinstance(child := content.get("children"), list):
            for child in content.get("children"):  # type: ignore
                child_either = ClasseqClade.from_dict(content=child)

                if child_either.is_left:
                    return child_either

                children.append(child_either.value)

        return right(
            cls(
                name=content.get("name"),  # type: ignore
                id=UUID(content.get("id")),
                type=eval(content.get("type")),  # type: ignore
                support=content.get("support"),
                branch_length=content.get("branch_length"),
                parent=(
                    UUID(parent) if (parent := content.get("parent")) else None
                ),
                children=tuple(children),
            )
        )

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def to_dict(self, omit_children: bool = False) -> dict[str, Any]:
        children = None

        if omit_children is False:
            children = (
                [child.to_dict() for child in self.children]
                if self.children is not None
                else None
            )

        return {
            "name": self.name,
            "id": self.id,
            "type": self.type,
            "support": self.support,
            "branch_length": self.branch_length,
            "parent": self.parent,
            "children": children,
        }

    def get_ingroup_clades(
        self,
    ) -> Either[c_exc.MappedErrors, list[Self]]:
        if (children := self.children) is None:
            return left(
                c_exc.ExecutionError(
                    "Unable to find outgroups. Specified clade has no "
                    + f"children: {self.get_pretty_clade()}",
                    exp=True,
                    logger=LOGGER,
                )
            )

        return right([clade for clade in children if clade.is_internal()])

    def get_outgroup_clade(self) -> Either[c_exc.MappedErrors, list[Self]]:
        if self.is_root() is False:
            return left(
                c_exc.ExecutionError(
                    "Could not get outgroup. The current clade is not the "
                    + "tree `root`.",
                    logger=LOGGER,
                )
            )

        if (children := self.children) is None:
            return left(
                c_exc.ExecutionError(
                    "Unable to find outgroups. Specified clade has no "
                    + f"children: {self.get_pretty_clade()}",
                    exp=True,
                    logger=LOGGER,
                )
            )

        return right([clade for clade in children if clade.is_outgroup()])

    def is_root(self) -> bool:
        return self.type == NodeType.ROOT

    def is_outgroup(self) -> bool:
        return self.type == NodeType.OUTGROUP

    def is_internal(self) -> bool:
        return self.type == NodeType.INTERNAL

    def is_terminal(self) -> bool:
        return self.type == NodeType.TERMINAL

    def get_pretty_clade(self) -> str:
        if self.is_terminal() or self.is_outgroup():
            return f"type({self.type.name}): {self.name}"

        nodes = 0 if self.children is None else len(self.children)
        name = self.name if self.name is not None else self.id
        return f"type({self.type.name}): [{nodes} nodes] {self.support} {name}"

    def get_pretty_tree(self) -> str:
        str_repr: list[str] = []
        backbone_color = "\033[90m"
        escape = "\033[0m"

        PIPE = f"{backbone_color}│{escape}"
        ELBOW = f"{backbone_color}└──{escape}"
        TEE = f"{backbone_color}├──{escape}"
        PIPE_PREFIX = f"{backbone_color}│   {escape}"
        SPACE_PREFIX = f"{backbone_color}    {escape}"

        def __build_node_representations(
            clade: ClasseqClade, prefix: str
        ) -> None:
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
                    colored_id = f"\033[94m{child.name or child.id}{escape}"

                    if entries_length - 1 == index:
                        indent_sep = SPACE_PREFIX
                        node_sep = ELBOW

                    str_repr.append(
                        prefix + node_sep + f"{colored_support} {colored_id}"
                    )

                    str_repr.append(prefix + indent_sep + PIPE)

                    __build_node_representations(
                        clade=child, prefix=prefix + indent_sep
                    )

        str_repr.append(PIPE)
        __build_node_representations(clade=self, prefix="")

        return "\n".join(str_repr)
