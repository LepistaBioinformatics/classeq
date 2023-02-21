import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER


def train_from_single_phylogeny(
    train_source: ReferenceSet,
) -> Either[bool, c_exc.MappedErrors]:
    try:
        # ? --------------------------------------------------------------------
        # ? Validate args
        # ? --------------------------------------------------------------------

        print(train_source)

        # ? --------------------------------------------------------------------
        # ? Recursively train for clades
        # ? --------------------------------------------------------------------

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(True)

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))
