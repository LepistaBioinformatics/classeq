from pathlib import Path
from typing import List

from Bio.Phylo.BaseTree import Tree

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER

from .collapse_low_supported_nodes import collapse_low_supported_nodes
from .parse_and_reroot_tree import parse_and_reroot_tree


def load_and_sanitize_phylogeny(
    source_file_path: Path,
    outgroups: List[str],
    support_value_cutoff: int = 99,
) -> Either[Tree, c_exc.MappedErrors]:
    """Load phylogenetic tree into memory.

    The loading process includes the literal loading using BioPython's Phylo
    module with further reroot and sanitizing. The sanitizing process include
    the collapsing of low supported branches.

    Args:
        source_file_path (Path): The system path of the input file.
        outgroups (List[str]): A list of outgroups to reroot the phylogenetic
            tree.
        support_value_cutoff (int, optional): A cutoff for collapsing internal
            nodes of the phylogenetic tree. Defaults to 99.

    Returns:
        Either[Tree, c_exc.MappedErrors]: A phylogenetic tree in the BioPython's
            format cas the use-case was successfully executed. A default
            MappedError otherwise.
    """

    try:
        # ? --------------------------------------------------------------------
        # ? Validate args
        # ? --------------------------------------------------------------------

        if not source_file_path.is_file():
            return left(
                c_exc.InvalidArgumentError(
                    f"Input file not exists: {source_file_path}"
                )
            )

        # ? --------------------------------------------------------------------
        # ? Load phylogenetic tree
        #
        # Load tree and extract out-group terminals.
        #
        # ? --------------------------------------------------------------------

        rooted_tree_either: Either = parse_and_reroot_tree(
            source_file_path, outgroups
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

        LOGGER.debug(f"rooted_tree: {rooted_tree}")
        LOGGER.debug(f"rooted_tree is rooted: {rooted_tree.rooted}")

        # ? --------------------------------------------------------------------
        # ? Sanitize tree
        #
        # Remove nodes with low phylogenetic support.
        #
        # ? --------------------------------------------------------------------

        sanitized_tree_either: Either = collapse_low_supported_nodes(
            rooted_tree, support_value_cutoff=support_value_cutoff
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

        LOGGER.debug(f"sanitized_tree: {sanitized_tree}")

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(sanitized_tree)

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))
