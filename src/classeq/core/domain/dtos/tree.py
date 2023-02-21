from pathlib import Path
from typing import List, Self

from attr import define, field
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.msa_source_format import MsaSourceFormatEnum
from classeq.core.domain.dtos.tree_source_format import TreeSourceFormatEnum
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER, TEMP_INPUT_FILE_SUFFIX


@define(kw_only=True)
class TreeSource:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    source_file_path: Path = field()
    tree_headers: List[str] = field()
    sanitized_tree: Tree | None = field(default=None)

    # ? ------------------------------------------------------------------------
    # ? Public class methods
    # ? ------------------------------------------------------------------------

    @classmethod
    def new(
        cls,
        source_file_path: Path,
        format: MsaSourceFormatEnum,
        outgroups: List[str],
        support_value_cutoff: int = 99,
    ) -> Either[Self, c_exc.MappedErrors]:
        try:
            if not source_file_path.is_file():
                return left(
                    c_exc.InvalidArgumentError(
                        f"Invalid path: {source_file_path}"
                    )
                )

            # ? ----------------------------------------------------------------
            # ? Load phylogenetic tree
            #
            # Load tree and extract out-group terminals.
            #
            # ? ----------------------------------------------------------------

            rooted_tree_either: Either = cls.__parse_and_reroot_tree(
                source_file_path,
                outgroups,
                format,
            )

            if rooted_tree_either.is_left:
                return left(
                    c_exc.UseCaseError(
                        "Unexpected error on parse phylogenetic tree.",
                        prev=rooted_tree_either.value,
                        logger=LOGGER,
                    )
                )

            rooted_tree: Tree = rooted_tree_either.value

            # ? ----------------------------------------------------------------
            # ? Sanitize tree
            #
            # Remove nodes with low phylogenetic support.
            #
            # ? ----------------------------------------------------------------

            sanitized_tree_either: Either = cls.__collapse_low_supported_nodes(
                rooted_tree,
                support_value_cutoff=support_value_cutoff,
            )

            if sanitized_tree_either.is_left:
                return left(
                    c_exc.UseCaseError(
                        "Unexpected error on sanitize phylogenetic tree.",
                        prev=sanitized_tree_either.value,
                        logger=LOGGER,
                    )
                )

            sanitized_tree: Tree = sanitized_tree_either.value

            # ? ----------------------------------------------------------------
            # ? Persist sanitized tree
            # ? ----------------------------------------------------------------

            cleaned_file_path = source_file_path.parent.joinpath(
                "".join(
                    [
                        source_file_path.stem,
                        ".",
                        TEMP_INPUT_FILE_SUFFIX,
                        source_file_path.suffix,
                    ]
                )
            )

            LOGGER.info(
                f"Sanitized TREE would be persisted to: {cleaned_file_path}"
            )

            Phylo.write(
                sanitized_tree,
                cleaned_file_path,
                format=format.value,
            )

            # ? ----------------------------------------------------------------
            # ? Return a positive response
            # ? ----------------------------------------------------------------

            return right(
                cls(
                    source_file_path=cleaned_file_path,
                    tree_headers=[
                        h.name for h in sanitized_tree.get_terminals()
                    ],
                    sanitized_tree=sanitized_tree,
                )
            )

        except Exception as exc:
            return left(c_exc.CreationError(exc, logger=LOGGER))

    # ? ------------------------------------------------------------------------
    # ? Private static methods
    # ? ------------------------------------------------------------------------

    @staticmethod
    def __parse_and_reroot_tree(
        source_file_path: Path,
        outgroups: List[str],
        format: TreeSourceFormatEnum,
    ) -> Either[Tree, c_exc.MappedErrors]:
        try:
            raw_tree: Tree = Phylo.read(source_file_path, format.value)

            terminals: List[Clade] = raw_tree.get_terminals()

            tree_outgroups = [
                clade for clade in terminals if clade.name in outgroups
            ]

            if not all([out.name in outgroups for out in tree_outgroups]):
                return left(
                    c_exc.UseCaseError(
                        f"Not all specified outgroups exists in phylogeny: {outgroups}",
                        exp=True,
                        logger=LOGGER,
                    )
                )

            raw_tree.root_with_outgroup([{"name": o} for o in outgroups])

            return right(raw_tree)

        except Exception as exc:
            return left(c_exc.UseCaseError(exc, logger=LOGGER))

    @staticmethod
    def __collapse_low_supported_nodes(
        rooted_tree: Tree,
        support_value_cutoff: int = 99,
    ) -> Either[Tree, c_exc.MappedErrors]:
        try:
            if rooted_tree.root is False:
                return left(
                    c_exc.InvalidArgumentError("Un-rooted trees not allowed.")
                )

            rooted_tree.collapse_all(
                lambda c: c.confidence is not None
                and c.confidence < support_value_cutoff
            )

            return right(rooted_tree)

        except Exception as exc:
            return left(c_exc.UseCaseError(exc, logger=LOGGER))
