from functools import reduce

import clean_base.exceptions as c_exc
from clean_base.either import Either, left, right

from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.priors import LabeledPriors
from classeq.settings import LOGGER

from ._calculate_probability_of_group_contains_kmer import (
    calculate_probability_of_group_contains_kmer,
)


def calculate_clade_adherence(
    labeled_priors: LabeledPriors,
    target_kmers: list[str],
    kmer_indices: KmersInverseIndices,
    total_length: int,
) -> Either[c_exc.MappedErrors, float]:
    """Calculate the probability of a sequence belongs to a clade.

    Description:
        This function calculates the probability of a sequence belongs to a
        clade. The probability is calculated by the joint probability of the
        sequence contains all the kmers of the clade.

    Args:
        labeled_priors (LabeledPriors): The labeled priors.
        target_kmers (list[str]): The target sequence kmers.
        kmer_indices (KmersInverseIndices): The kmer indices.
        total_length (int): The total length of the sequences.

    Returns:
        Either[c_exc.MappedErrors, float]: The probability of the sequence
            belongs to the clade.

    Raises:
        c_exc.UseCaseError: If the kmer is not indexed.
    """
    try:
        joint_probability_units: list[float] = []

        for kmer in target_kmers:
            if (prior := labeled_priors.priors.get(kmer)) is None:
                # LOGGER.debug(">>>>>> CONTINUE <<<<<<<<<<<")
                continue

            if (kmer_index := kmer_indices.index_of(kmer)) is None:
                return left(
                    c_exc.UseCaseError(
                        "Unexpected error on calculate outgroup joint "
                        + f"probability. `{kmer}` is not indexed.",
                        logger=LOGGER,
                    )
                )

            current_index = [
                i
                for i in labeled_priors.labels
                if kmer_indices.indices[kmer_index].contains(i)
            ]

            joint_probability_units.append(
                calculate_probability_of_group_contains_kmer(
                    prior=prior,
                    sequences_with_kmer=len(current_index),
                    total_sequences=total_length,
                )
            )

        return right(reduce(lambda i, j: i * j, joint_probability_units))

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))
