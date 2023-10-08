import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.ordered_tuple import OrderedTuple
from classeq.settings import LOGGER


def fetch_kmers_indices(
    kmer_indices: KmersInverseIndices,
    sequence_codes: list[int],
) -> Either[c_exc.MappedErrors, set[str]]:
    """Fetch kmer indices for a given set of sequence codes.

    Args:
        kmer_indices (KmersInverseIndices): Kmer indices of the MSA.
        sequence_codes (list[int]): Sequence codes to be used as filter.

    Returns:
        Either[c_exc.MappedErrors, set[str]]: Either a set of kmer indices
            or a `classeq.core.domain.utils.exceptions.MappedErrors`
            instance.

    Raises:
        c_exc.UseCaseError: If any error occurred.

    """

    try:
        return right(
            OrderedTuple(
                {
                    index.kmer
                    for index in kmer_indices.indices
                    if not set(index.records).isdisjoint(sequence_codes)
                }
            )
        )

    except Exception as exc:
        raise c_exc.UseCaseError(exc, logger=LOGGER)()
