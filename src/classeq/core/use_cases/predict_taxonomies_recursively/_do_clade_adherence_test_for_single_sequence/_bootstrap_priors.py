from random import choices
from typing import Any

import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.ordered_tuple import OrderedTuple
from classeq.core.domain.dtos.priors import NoiseGroupLabeledPriors
from classeq.settings import LOGGER

from .._calculate_clade_adherence_with_bootstrap import (
    AdherenceStatus,
    calculate_clade_adherence_test,
)


def bootstrap_priors(
    query_kmers: set[str],
    kmer_indices: KmersInverseIndices,
    bootstrap_replicates: int,
    sample_size: int,
    noise_group_priors: dict[str, float],
    noise_group_labels: set[int],
    **kwargs: Any,
) -> Either[c_exc.MappedErrors, tuple[list[float], list[int]]]:
    """Calculate bootstrap for adherence test.

    Args:
        query_kmers (set[str]): The query kmers.
        kmer_indices (KmersInverseIndices): The kmer indices.
        bootstrap_replicates (int): The number of bootstrap replicates.
        sample_size (int): The sample size.
        noise_group_priors (dict[str, float]): The noise group priors.
        noise_group_labels (set[str]): The noise group labels.

    Returns:
        Either[c_exc.MappedErrors, tuple[list[float], list[int]]]: The
            bootstrap results.

    Raises:
        c_exc.UseCaseError: If any of the arguments is not a list.
        c_exc.UseCaseError: If any of the sub-items of the arguments is
            not a list.

    """

    joint_probabilities: list[float] = []
    matching_kmers: list[int] = []

    try:
        for _ in range(bootstrap_replicates):
            priors_sample = choices(
                list(noise_group_priors.items()),
                k=sample_size,
            )

            noise_group = NoiseGroupLabeledPriors(
                labels=OrderedTuple(noise_group_labels),
                priors={item[0]: item[1] for item in priors_sample},
            )

            if (
                adherence_either := calculate_clade_adherence_test(
                    labeled_priors=noise_group,
                    query_kmers=query_kmers,
                    kmer_indices=kmer_indices,
                    total_length=len(noise_group.labels),
                    **kwargs,
                )
            ).is_left:
                return adherence_either

            adherence = adherence_either.value
            if adherence.status == AdherenceStatus.SUCCESS:
                joint_probabilities.append(adherence.joint_probability)
                matching_kmers.append(adherence.match_kmers)

        return right((joint_probabilities, matching_kmers))

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
