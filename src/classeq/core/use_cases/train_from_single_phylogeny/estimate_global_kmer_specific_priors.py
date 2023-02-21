from collections import defaultdict
from typing import DefaultDict

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.kmer_inverse_index import KmerIndex
from classeq.core.domain.dtos.msa import MsaSource
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER


def estimate_global_kmer_specific_priors(
    msa: MsaSource,
) -> Either[DefaultDict[str, float], c_exc.MappedErrors]:
    try:
        # ? --------------------------------------------------------------------
        # ? Validate input args
        # ? --------------------------------------------------------------------

        if not isinstance(msa, MsaSource):
            return left(
                c_exc.InvalidArgumentError(
                    "`msa` is not a instance of `MsaSource`.",
                    exp=True,
                    logger=LOGGER,
                )
            )

        # ? --------------------------------------------------------------------
        # ? Calculate by-kmer indices
        # ? --------------------------------------------------------------------

        output: DefaultDict[str, float] = defaultdict()

        for index in msa.kmers_indices.indices:
            index_response_either = __calculate_single_kmer_likelihood(
                total_records=len(msa.sequence_headers),
                kmers_index=index,
            )

            if index_response_either.is_left:
                return left(
                    c_exc.ExecutionError(
                        f"Unexpected error on calculate likelihood of kmer {index.kmer}.",
                        prev=index_response_either.value,
                        logger=LOGGER,
                    )
                )

            output[index.kmer] = index_response_either.value

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(output)

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))


def __calculate_single_kmer_likelihood(
    total_records: int,
    kmers_index: KmerIndex,
) -> Either[float, c_exc.MappedErrors]:
    try:
        return right((len(kmers_index.records) + 0.5) / (total_records + 1))

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))
