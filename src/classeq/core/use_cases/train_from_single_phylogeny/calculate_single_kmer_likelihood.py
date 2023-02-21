import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.kmer_inverse_index import KmerIndex
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER


def calculate_single_kmer_likelihood(
    total_records: int,
    kmers_index: KmerIndex,
) -> Either[float, c_exc.MappedErrors]:
    try:
        return right((len(kmers_index.records) + 0.5) / (total_records + 1))

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))
