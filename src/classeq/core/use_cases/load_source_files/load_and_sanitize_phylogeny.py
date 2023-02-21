from pathlib import Path
from typing import List

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.tree import TreeSource
from classeq.core.domain.dtos.tree_source_format import TreeSourceFormatEnum
from classeq.core.domain.utils.either import Either, left
from classeq.settings import LOGGER


def load_and_sanitize_phylogeny(
    source_file_path: Path,
    outgroups: List[str],
    format: TreeSourceFormatEnum,
    support_value_cutoff: int = 99,
) -> Either[TreeSource, c_exc.MappedErrors]:
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
        return TreeSource.new(
            source_file_path=source_file_path,
            format=format,
            outgroups=outgroups,
            support_value_cutoff=support_value_cutoff,
        )

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))
