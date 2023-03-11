import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.clade import CladeWrapper
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER


def perform_phylogenetic_adherence_test(
    target_sequence: str,
    seed_clade: CladeWrapper,
) -> Either[c_exc.MappedErrors, bool]:
    try:
        # ? --------------------------------------------------------------------
        # ? Validate args
        # ? --------------------------------------------------------------------

        # TODO: do implement

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(True)

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))
