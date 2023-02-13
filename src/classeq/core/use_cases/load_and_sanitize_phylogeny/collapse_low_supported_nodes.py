from Bio.Phylo.BaseTree import Tree

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER


def collapse_low_supported_nodes(
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
