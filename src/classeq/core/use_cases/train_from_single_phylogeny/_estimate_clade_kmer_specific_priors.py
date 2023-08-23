from collections import defaultdict

import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.domain.dtos.kmer_inverse_index import (
    KmerIndex,
    KmersInverseIndices,
)
from classeq.settings import LOGGER


def estimate_kmer_specific_priors_of_clade(
    kmer_indices: KmersInverseIndices,
    sequence_codes: list[int],
    corpus_size: int,
) -> Either[c_exc.MappedErrors, defaultdict[str, float]]:
    """Estimate clade specific priors for each kmer.

    Args:
        kmer_indices (KmersInverseIndices): Kmers inverse indices.
        sequence_codes (list[int]): Sequence codes.
        corpus_size (int): Corpus size.

    Returns:
        Either[c_exc.MappedErrors, DefaultDict[str, float]]: Either a
            UseCaseError or a dictionary with the kmer specific priors.

    Raises:
        UseCaseError: If the kmer index does not contains any sequence code.

    """

    try:
        kmers_priors_for_clade: defaultdict[str, float] = defaultdict()
        target_indices: set[KmerIndex] = set()

        for index in kmer_indices.indices:
            if set(index.records).isdisjoint(sequence_codes):
                continue

            target_indices.add(index)

        for index in iter(target_indices):
            n = len(set(index.records).intersection(sequence_codes))

            if n == 0:
                raise

            prior = (n + 0.5) / (corpus_size + 1)
            kmers_priors_for_clade[index.kmer] = prior

        return right(kmers_priors_for_clade)

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
