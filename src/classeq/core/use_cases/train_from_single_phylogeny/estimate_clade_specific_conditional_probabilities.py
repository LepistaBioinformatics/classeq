import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.msa import MsaSource
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER


def estimate_clade_specific_conditional_probabilities(
    msa: MsaSource,
) -> Either[bool, c_exc.MappedErrors]:
    try:
        # ? --------------------------------------------------------------------
        # ? Validate entries
        # ? --------------------------------------------------------------------

        if msa.kmers_indices is None:
            LOGGER.warning("Kmer probabilities not initialized. Doing it!")

            init_either = msa.initialize_kmer_indices()

            if init_either.is_left:
                return left(
                    c_exc.InvalidArgumentError(
                        "Unexpected error on initialize kmers prior probabilities.",
                        prev=init_either.value,
                        logger=LOGGER,
                    )
                )

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(True)

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))
