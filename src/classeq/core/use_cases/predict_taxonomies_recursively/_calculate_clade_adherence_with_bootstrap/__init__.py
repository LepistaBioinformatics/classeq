from functools import reduce

import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.priors import LabeledPriors
from classeq.settings import LOGGER

from ._calculate_probability_of_group_contains_kmer import (
    calculate_probability_of_group_contains_kmer,
)
from ._dtos import AdherenceResult, AdherenceStatus, AdherenceTestStrategy


def calculate_clade_adherence_test(
    labeled_priors: LabeledPriors,
    query_kmers: set[str],
    kmer_indices: KmersInverseIndices,
    total_length: int,
    adherence_strategy: AdherenceTestStrategy,
) -> Either[c_exc.MappedErrors, AdherenceResult]:
    """Calculate the probability of a sequence belongs to a clade.

    Description:
        This function calculates the probability of a sequence belongs to a
        clade. The probability is calculated by the joint probability of the
        sequence contains all the kmers of the clade.

    Args:
        labeled_priors (LabeledPriors): The labeled priors. target_kmers
        (list[str]): The target sequence kmers. kmer_indices
        (KmersInverseIndices): The kmer indices. total_length (int): The total
        length of the sequences.

    Returns:
        Either[c_exc.MappedErrors, AdherenceResult]: The probability of the
            sequence belongs to the clade wrapped into a AdherenceResult object.

    Raises:
        c_exc.UseCaseError: If the kmer is not indexed.

    """

    try:
        full_joint_probability_units: list[float] = []
        labeled_priors_keys = labeled_priors.priors.keys()
        kmers_intersection = query_kmers.intersection(set(labeled_priors_keys))

        if len(kmers_intersection) == 0:
            return right(
                AdherenceResult(
                    match_kmers=0,
                    query_kmers_size=len(query_kmers),
                    subject_kmers_size=len(labeled_priors_keys),
                    status=AdherenceStatus.NOT_ENOUGH_PRIORS,
                    strategy=adherence_strategy,
                )
            )

        if adherence_strategy == AdherenceTestStrategy.KMERS_INTERSECTION:
            return right(
                AdherenceResult(
                    match_kmers=len(kmers_intersection),
                    query_kmers_size=len(query_kmers),
                    subject_kmers_size=len(labeled_priors_keys),
                    status=AdherenceStatus.SUCCESS,
                    strategy=adherence_strategy,
                )
            )

        for kmer in kmers_intersection:
            if (prior := labeled_priors.priors.get(kmer)) is None:
                return c_exc.UseCaseError(
                    "Unexpected error on calculate joint probability of "
                    + f"`{kmer}`. There is not indexed.",
                    logger=LOGGER,
                )()

            if (kmer_index := kmer_indices.index_of(kmer)) is None:
                return c_exc.UseCaseError(
                    "Unexpected error on calculate outgroup joint "
                    + f"probability. `{kmer}` is not indexed.",
                    logger=LOGGER,
                )()

            current_index = [
                label
                for label in labeled_priors.labels
                if kmer_indices.indices[kmer_index].contains(label)
            ]

            full_joint_probability_units.append(
                calculate_probability_of_group_contains_kmer(
                    prior=prior,
                    sequences_with_kmer=len(current_index),
                    total_sequences=total_length,
                )
            )

        full_joint_probability = reduce(
            lambda i, j: i * j, full_joint_probability_units
        )

        return right(
            AdherenceResult(
                match_kmers=len(kmers_intersection),
                query_kmers_size=len(query_kmers),
                subject_kmers_size=len(labeled_priors_keys),
                status=AdherenceStatus.SUCCESS,
                joint_probability=full_joint_probability,
            )
        )

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
