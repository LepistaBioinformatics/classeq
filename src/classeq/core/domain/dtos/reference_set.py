from typing import DefaultDict, Set
from uuid import UUID, uuid4

from attr import field, frozen
from Bio.Phylo.BaseTree import Clade, Tree

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.clade import CladeWrapper, NodeType
from classeq.core.domain.dtos.msa import MsaSource
from classeq.core.domain.dtos.tree import TreeSource
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER


@frozen(kw_only=True)
class ReferenceSet:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    tree: TreeSource = field()
    msa: MsaSource = field()
    labels_map: DefaultDict[str, int] = field()

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def get_linear_tree(self) -> Either[Set[CladeWrapper], c_exc.MappedErrors]:
        """Convert the original Bio.Phylo.BaseTree.Tree to a linear version
        composed of a set of `CladeWrapper` elements.
        """

        def __collect_clades_recursively(
            parent: UUID | None,
            clade: Clade,
        ) -> None:
            """Recursively collect nodes an parse to CladeWrappers.

            Args:
                parent (UUID | None): The UUID of the parent node.
                clade (Clade): The target clade.
            """

            current_clade_identifier = uuid4()

            linear_tree.add(
                CladeWrapper(
                    id=current_clade_identifier,
                    name=clade.name,
                    type=(
                        NodeType.TERMINAL
                        if clade.is_terminal()
                        else NodeType.INTERNAL
                    ),
                    parent=parent,
                    support=clade.confidence,
                )
            )

            if len(clade.clades) > 0:
                for clade in clade.clades:
                    __collect_clades_recursively(
                        parent=current_clade_identifier,
                        clade=clade,
                    )

        try:
            tree: Tree = self.tree.sanitized_tree
            root: Clade = tree.root

            linear_tree: Set[CladeWrapper] = set()
            linear_tree.add(
                CladeWrapper(
                    name=NodeType.ROOT.value,
                    type=NodeType.ROOT,
                )
            )

            __collect_clades_recursively(
                parent=None,
                clade=root,
            )

            return right(Set[CladeWrapper])

        except Exception as exc:
            return left(c_exc.UseCaseError(exc, logger=LOGGER))
