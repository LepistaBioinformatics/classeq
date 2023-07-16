from functools import reduce

import clean_base.exceptions as c_exc
from clean_base.either import Either, left, right

from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.priors import (
    IngroupCladePriors,
    LabeledPriors,
    OutgroupCladePriors,
    PriorGroup,
)
from classeq.settings import DEFAULT_KMER_SIZE, LOGGER


def do_clade_adherence_test_for_single_sequence(
    target_sequence: str,
    clade_priors: IngroupCladePriors | OutgroupCladePriors,
    kmer_indices: KmersInverseIndices,
    total_length: int,
) -> Either[c_exc.MappedErrors, dict[PriorGroup, float]]:
    """Calculate the probability of a sequence belongs to a clade.

    Description:
        This function calculates the probability of a sequence belongs to a
        clade. The probability is calculated by the joint probability of the
        sequence contains all the kmers of the clade.

    Args:
        target_sequence (str): The sequence to be tested.
        clade_priors (IngroupCladePriors | OutgroupCladePriors): The clade
            priors.
        kmer_indices (KmersInverseIndices): The kmer indices.
        total_length (int): The total length of the sequences.

    Returns:
        Either[c_exc.MappedErrors, dict[PriorGroup, float]]: The probability of
            the sequence belongs to the clade.

    Raises:
        c_exc.InvalidArgumentError: If any of the arguments is not a list.
        c_exc.InvalidArgumentError: If any of the sub-items of the arguments is
            not a list.

    """
    try:
        # ? --------------------------------------------------------------------
        # ? Validate entries
        # ? --------------------------------------------------------------------

        if not any(
            [
                isinstance(clade_priors, IngroupCladePriors),
                isinstance(clade_priors, OutgroupCladePriors),
            ]
        ):
            return left(
                c_exc.InvalidArgumentError(
                    f"argument `{clade_priors}` should be any of "
                    + f"`{IngroupCladePriors}` or `{OutgroupCladePriors}`.",
                    exp=True,
                    logger=LOGGER,
                )
            )

        for observed, desired in [
            (target_sequence, str),
            (kmer_indices, KmersInverseIndices),
        ]:
            if not isinstance(observed, desired):
                return left(
                    c_exc.InvalidArgumentError(
                        f"argument `{observed}` should be a instance of `{desired}`.",
                        exp=True,
                        logger=LOGGER,
                    )
                )

        # ? --------------------------------------------------------------------
        # ? Get sequence kmers
        # ? --------------------------------------------------------------------

        target_sequence_kmers = [
            kmer
            for kmer in KmersInverseIndices.generate_kmers(
                dna_sequence=target_sequence.upper(),
                k_size=DEFAULT_KMER_SIZE,
            )
        ]

        # ? --------------------------------------------------------------------
        # ? Calculate joint probability units
        # ? --------------------------------------------------------------------

        clade_adherence_stats: dict[PriorGroup, float] = {}

        if isinstance(clade_priors, IngroupCladePriors):
            for group in clade_priors.priors:
                adherence_either = __calculate_clade_adherence(
                    labeled_priors=group,
                    target_kmers=target_sequence_kmers,
                    kmer_indices=kmer_indices,
                    # total_length=sum(
                    #     len(i.labels) for i in clade_priors.priors
                    # ),
                    total_length=total_length,
                )

                if adherence_either.is_left:
                    return adherence_either

                clade_adherence_stats.update(
                    {group.group: adherence_either.value}
                )

        if isinstance(clade_priors, OutgroupCladePriors):
            adherence_either = __calculate_clade_adherence(
                labeled_priors=clade_priors.priors,
                target_kmers=target_sequence_kmers,
                kmer_indices=kmer_indices,
                total_length=len(clade_priors.priors.labels),
            )

            if adherence_either.is_left:
                return adherence_either

            clade_adherence_stats.update(
                {clade_priors.priors.group: adherence_either.value}
            )

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(clade_adherence_stats)

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))


def __calculate_clade_adherence(
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
                __calculate_probability_of_group_contains_kmer(
                    prior=prior,
                    sequences_with_kmer=len(current_index),
                    # total_sequences=len(labeled_priors.labels),
                    total_sequences=total_length,
                )
            )

        return right(reduce(lambda i, j: i * j, joint_probability_units))

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))


def __calculate_probability_of_group_contains_kmer(
    prior: float,
    sequences_with_kmer: int,
    total_sequences: int,
) -> float:
    """Calculate the probability of a group contains a kmer.

    Description:
        This function calculates the probability of a group contains a kmer.
        The probability is calculated by the following formula:
        P(wi|G) = [ m(wi) + Pi ] / ( M + 1 ) where m(wi) is the number of the
        group sequences containing word wi, Pi is the prior probability for the
        kmer in the overall training dataset, and M is the dataset size for the
        target training group.

    Args:
        prior (float): The prior probability for the kmer in the overall
            training dataset.
        sequences_with_kmer (int): The number of the group sequences
            containing word wi.
        total_sequences (int): The dataset size for the target training group.

    Returns:
        float: The probability of a group contains a kmer.

    Raises:
        Exception: If `total_sequences` is greater than `sequences_with_kmer`.

    """

    if sequences_with_kmer > total_sequences:
        raise Exception(
            "`total_sequences` could not be greater than `sequences_with_kmer`"
        )

    # ? ------------------------------------------------------------------------
    # ? Calculate
    #
    # The original formulation for the probability that a group contains a kmer
    # is:
    #
    # P(wi|G) = [ m(wi) + Pi ] / ( M + 1 )
    #
    # m(wi) = be the number of the group sequences containing word wi (argument
    # `sequences_with_kmer` of these function).
    #
    # Pi = The prior probability for the kmer in the overall training dataset
    # (argument `prior` of these function).
    #
    # M = The dataset size for the target training group (argument
    # `total_sequences` of these function).
    #
    # ? ------------------------------------------------------------------------

    return (sequences_with_kmer + prior) / (total_sequences + 1)
