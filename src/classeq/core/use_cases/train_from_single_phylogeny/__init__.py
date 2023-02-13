from Bio.Phylo.BaseTree import Tree

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER


def train_from_single_phylogeny(tree: Tree) -> Either[bool, c_exc.MappedErrors]:
    try:
        # ? --------------------------------------------------------------------
        # ? Validate args
        # ? --------------------------------------------------------------------

        # ? --------------------------------------------------------------------
        # ? Recursively train for clades
        # ? --------------------------------------------------------------------

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(True)

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))
