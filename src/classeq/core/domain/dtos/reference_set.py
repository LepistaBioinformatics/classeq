from collections import defaultdict
from typing import Any, DefaultDict, Self
from uuid import UUID, uuid4

from attr import field, define
from Bio.Phylo.BaseTree import Clade, Tree

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.clade import CladeWrapper, NodeType
from classeq.core.domain.dtos.msa import MsaSource
from classeq.core.domain.dtos.tree import TreeSource
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER


@define(kw_only=True)
class ReferenceSet:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    tree: TreeSource = field()
    msa: MsaSource = field()
    labels_map: DefaultDict[str, int] = field()

    # ? ------------------------------------------------------------------------
    # ? Public class methods
    # ? ------------------------------------------------------------------------

    @classmethod
    def from_dict(
        cls,
        content: dict[str, Any],
    ) -> Either[c_exc.MappedErrors, Self]:
        for key in [
            "tree",
            "msa",
            "labels_map",
        ]:
            if key not in content:
                return left(
                    c_exc.InvalidArgumentError(
                        f"Invalid content detected on parse `{ReferenceSet}`. "
                        f"{key}` key is empty.",
                        logger=LOGGER,
                    )
                )

        tree_either = TreeSource.from_dict(content=content.get("tree"))

        if tree_either.is_left:
            return tree_either

        msa_either = MsaSource.from_dict(content=content.get("msa"))

        if msa_either.is_left:
            return msa_either

        return right(
            cls(
                tree=tree_either.value,
                msa=msa_either.value,
                labels_map=defaultdict(int, content.get("labels_map")),  # type: ignore
            )
        )

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def get_hierarchical_tree(
        self,
    ) -> Either[c_exc.MappedErrors, CladeWrapper]:
        try:
            # ? ----------------------------------------------------------------
            # ? Generate the linear tree
            # ? ----------------------------------------------------------------

            linear_tree_either = self.get_linear_tree()

            if linear_tree_either.is_left:
                return linear_tree_either

            linear_tree: set[CladeWrapper] = linear_tree_either.value

            if len([i for i in linear_tree if i.is_root()]) != 1:
                return left(
                    c_exc.ExecutionError(
                        "More than one root node found in linear tree.",
                        logger=LOGGER,
                    )
                )

            # ? ----------------------------------------------------------------
            # ? Get the tree seed
            # ? ----------------------------------------------------------------

            try:
                root = next(i for i in iter(linear_tree) if i.is_root())
            except StopIteration:
                return left(
                    c_exc.ExecutionError(
                        "Root node not present in linear tree.",
                        logger=LOGGER,
                    )
                )

            # ? ----------------------------------------------------------------
            # ? Start the tree expansion
            # ? ----------------------------------------------------------------

            seed_tree: CladeWrapper = root

            def __expand_tree(clade: CladeWrapper) -> None:
                for child in [i for i in linear_tree if i.parent == clade.id]:
                    if clade.children is None:
                        clade.children = tuple()

                    clade.children = clade.children + (child,)

                    if child.is_internal() or child.is_outgroup():
                        __expand_tree(clade=child)

            __expand_tree(clade=seed_tree)

            # ? ----------------------------------------------------------------
            # ? Return a positive response
            # ? ----------------------------------------------------------------

            return right(seed_tree)

        except Exception as exc:
            return left(c_exc.ExecutionError(exc, logger=LOGGER))

    def get_linear_tree(
        self,
    ) -> Either[c_exc.MappedErrors, set[CladeWrapper]]:
        """Convert the original Bio.Phylo.BaseTree.Tree to a linear version
        composed of a set of `CladeWrapper` elements.
        """

        def __determine_node_type(clade: Clade) -> NodeType:
            if clade.is_terminal():
                if clade.name in self.tree.outgroups:
                    return NodeType.OUTGROUP
                return NodeType.TERMINAL
            return NodeType.INTERNAL

        def __collect_clades_recursively(
            parent: UUID | None,
            clade: Clade,
            is_root: bool = False,
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
                    name=(
                        NodeType.ROOT.value if is_root is True else clade.name
                    ),
                    type=(
                        NodeType.ROOT
                        if is_root is True
                        else __determine_node_type(clade=clade)
                    ),
                    parent=parent,
                    support=clade.confidence,
                    branch_length=clade.branch_length,
                )
            )

            if len(clade.clades) > 0:
                for clade in clade.clades:
                    __collect_clades_recursively(
                        parent=current_clade_identifier,
                        clade=clade,
                    )

        try:
            if self.tree.sanitized_tree is None:
                tree_either = self.tree.parse_and_reroot_tree(
                    source_file_path=self.tree.source_file_path,
                    outgroups=self.tree.outgroups,
                    format=self.tree.tree_format,
                )

                if tree_either.is_left:
                    return tree_either

                self.tree.sanitized_tree = tree_either.value

            tree: Tree = self.tree.sanitized_tree
            root: Clade = tree.root

            linear_tree: set[CladeWrapper] = set()

            __collect_clades_recursively(
                parent=None,
                clade=root,
                is_root=True,
            )

            return right(linear_tree)

        except Exception as exc:
            return left(c_exc.ExecutionError(exc, logger=LOGGER))
