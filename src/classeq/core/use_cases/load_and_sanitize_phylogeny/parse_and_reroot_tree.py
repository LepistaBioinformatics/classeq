from pathlib import Path
from typing import List

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER


def parse_and_reroot_tree(
    source_file_path: Path,
    outgroups: List[str],
) -> Either[Tree, c_exc.MappedErrors]:
    try:
        raw_tree: Tree = Phylo.read(source_file_path, "newick")

        terminals: List[Clade] = raw_tree.get_terminals()

        tree_outgroups = [
            clade for clade in terminals if clade.name in outgroups
        ]

        if not all([out.name in outgroups for out in tree_outgroups]):
            return left(c_exc.UseCaseError("", logger=LOGGER))

        raw_tree.root_with_outgroup([{"name": o for o in outgroups}])

        return right(raw_tree)

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))
