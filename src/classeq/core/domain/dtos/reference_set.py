from collections import defaultdict
from copy import deepcopy
from typing import Any, DefaultDict, Self
from uuid import UUID

import clean_base.exceptions as c_exc
from attr import define, field
from clean_base.either import Either, right

from classeq.core.domain.dtos.biopython_wrappers import (
    ExtendedBioPythonClade,
    ExtendedBioPythonTree,
)
from classeq.core.domain.dtos.clade import ClasseqClade, NodeType
from classeq.core.domain.dtos.msa import MsaSource
from classeq.core.domain.dtos.tree import ClasseqTree
from classeq.settings import LOGGER


@define(kw_only=True)
class ReferenceSet:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    tree: ClasseqTree = field()
    msa: MsaSource = field()
    labels_map: DefaultDict[str, int] = field()
    linear_tree: tuple[ClasseqClade, ...] | None = field(default=None)

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
            "linear_tree",
        ]:
            if key not in content:
                return c_exc.DadaTransferObjectError(
                    f"Invalid content detected on parse `{ReferenceSet}`. "
                    f"{key}` key is empty.",
                    logger=LOGGER,
                )()

        if (
            tree_either := ClasseqTree.from_dict(content=content.pop("tree"))
        ).is_left:
            return tree_either

        if (
            msa_either := MsaSource.from_dict(content=content.pop("msa"))
        ).is_left:
            return msa_either

        linear_tree: list[ClasseqClade] | None = None
        if isinstance(linear_tree_content := content.pop("linear_tree"), list):
            for linear_tree_unit in linear_tree_content:
                if (
                    unit_either := ClasseqClade.from_dict(linear_tree_unit)
                ).is_left:
                    return unit_either

                if linear_tree is None:
                    linear_tree = []

                linear_tree.append(unit_either.value)

        return right(
            cls(
                tree=tree_either.value,
                msa=msa_either.value,
                labels_map=defaultdict(int, content.pop("labels_map")),  # type: ignore
                linear_tree=(
                    linear_tree if linear_tree is None else tuple(linear_tree)
                ),
            )
        )

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def get_hierarchical_tree(
        self,
        force_reload_tree: bool = False,
    ) -> Either[c_exc.MappedErrors, ClasseqClade]:
        try:
            # ? ----------------------------------------------------------------
            # ? Generate the linear tree
            # ? ----------------------------------------------------------------

            if self.linear_tree is None or force_reload_tree is True:
                linear_tree_either = self.build_linear_tree(force_reload_tree)

                if linear_tree_either.is_left:
                    return linear_tree_either

            # The `linear_tree` attribute should exists after execution of the
            # `build_linear_tree` method. Don't remove these check.
            if (linear_tree := self.linear_tree) is None:
                return c_exc.DadaTransferObjectError(
                    "`build_linear_tree` method is maybe not working. "
                    + f"Attribute `linear_tree` of {Self} is `None` after "
                    + "execute it. The expected behavior is to be of type"
                    + f"`{tuple[ClasseqClade, ...]}`.",
                    logger=LOGGER,
                )()

            if len([i for i in linear_tree if i.is_root()]) != 1:
                return c_exc.DadaTransferObjectError(
                    "More than one root node found in linear tree.",
                    logger=LOGGER,
                )()

            # ? ----------------------------------------------------------------
            # ? Get the tree seed
            # ? ----------------------------------------------------------------

            try:
                root = next(i for i in iter(linear_tree) if i.is_root())
            except StopIteration:
                return c_exc.DadaTransferObjectError(
                    "Root node not present in linear tree.",
                    logger=LOGGER,
                )()

            # ? ----------------------------------------------------------------
            # ? Start the tree expansion
            # ? ----------------------------------------------------------------

            seed_tree: ClasseqClade = deepcopy(root)

            def __expand_tree(
                clade: ClasseqClade,
                linear_tree: tuple[ClasseqClade, ...],
            ) -> None:
                for child in [i for i in linear_tree if i.parent == clade.id]:
                    if clade.children is None:
                        clade.children = tuple()

                    clade.children = clade.children + (child,)

                    if child.is_internal():
                        __expand_tree(clade=child, linear_tree=linear_tree)

            __expand_tree(clade=seed_tree, linear_tree=linear_tree)

            # ? ----------------------------------------------------------------
            # ? Return a positive response
            # ? ----------------------------------------------------------------

            return right(seed_tree)

        except Exception as exc:
            return c_exc.DadaTransferObjectError(exc, logger=LOGGER)()

    def build_linear_tree(
        self,
        force_reload_tree: bool = False,
    ) -> Either[c_exc.MappedErrors, bool]:
        """Convert the original Bio.Phylo.BaseTree.Tree to a linear version
        composed of a set of `CladeWrapper` elements.
        """

        def __discover_node_type(clade: ExtendedBioPythonClade) -> NodeType:
            if clade.is_terminal():
                if clade.name in self.tree.outgroups:
                    return NodeType.OUTGROUP
                return NodeType.TERMINAL
            return NodeType.INTERNAL

        def __collect_clades_recursively(
            parent: UUID | None,
            clade: ExtendedBioPythonClade,
            is_root: bool = False,
        ) -> None:
            """Recursively collect nodes an parse to CladeWrappers.

            Args:
                parent (UUID | None): The UUID of the parent node.
                clade (Clade): The target clade.
            """

            linear_tree.add(
                ClasseqClade(
                    id=clade._id,
                    name=(
                        NodeType.ROOT.value
                        if is_root is True
                        else self.labels_map.get(clade.name)
                    ),
                    type=(
                        NodeType.ROOT
                        if is_root is True
                        else __discover_node_type(clade=clade)
                    ),
                    parent=parent,
                    support=clade.confidence,
                    branch_length=clade.branch_length,
                )
            )

            if len(clade.clades) > 0:
                for children in clade.clades:
                    __collect_clades_recursively(
                        parent=clade._id,
                        clade=children,
                    )

        try:
            if self.tree.sanitized_tree is None or force_reload_tree is True:
                if (
                    tree_either := self.tree.parse_and_reroot_phylojson_tree(
                        phylojson_file_path=self.tree.phylojson_file_path,
                        outgroups=self.tree.outgroups,
                    )
                ).is_left:
                    return tree_either

                self.tree.sanitized_tree = tree_either.value

            tree: ExtendedBioPythonTree = deepcopy(self.tree.sanitized_tree)
            root: ExtendedBioPythonClade = deepcopy(tree.root)

            linear_tree: set[ClasseqClade] = set()

            __collect_clades_recursively(
                parent=None,
                clade=root,
                is_root=True,
            )

            self.linear_tree = tuple(linear_tree)

            return right(True)

        except Exception as exc:
            return c_exc.DadaTransferObjectError(exc, logger=LOGGER)()
