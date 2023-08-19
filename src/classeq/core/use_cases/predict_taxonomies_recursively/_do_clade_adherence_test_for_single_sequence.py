from copy import copy
from random import choices, seed
from typing import Any

import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.priors import (
    IngroupCladePriors,
    IngroupLabeledPriors,
    LabeledPriors,
    NoiseGroupLabeledPriors,
    OutgroupCladePriors,
    PriorGroup,
    SisterGroupLabeledPriors,
)
from classeq.settings import LOGGER

from ._calculate_clade_adherence_with_bootstrap import (
    AdherenceResult,
    AdherenceStatus,
    AdherenceTestStrategy,
    calculate_clade_adherence_test,
)


def do_clade_adherence_test_for_single_sequence(
    query_kmers: set[str],
    clade_priors: IngroupCladePriors | OutgroupCladePriors,
    kmer_indices: KmersInverseIndices,
    calculate_bootstrap: bool = False,
    bootstrap_replicates: int = 100,
    adherence_strategy: AdherenceTestStrategy = AdherenceTestStrategy(None),
) -> Either[c_exc.MappedErrors, dict[PriorGroup, AdherenceResult]]:
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
        c_exc.UseCaseError: If any of the arguments is not a list.
        c_exc.UseCaseError: If any of the sub-items of the arguments is
            not a list.

    """

    try:
        seed(999_999_999)

        # ? --------------------------------------------------------------------
        # ? Validate entries
        # ? --------------------------------------------------------------------

        if not any(
            [
                isinstance(clade_priors, IngroupCladePriors),
                isinstance(clade_priors, OutgroupCladePriors),
            ]
        ):
            return c_exc.UseCaseError(
                f"argument `{clade_priors}` should be any of "
                + f"`{IngroupCladePriors}` or `{OutgroupCladePriors}`.",
                exp=True,
                logger=LOGGER,
            )()

        for observed, desired in [
            (query_kmers, set),
            (kmer_indices, KmersInverseIndices),
        ]:
            if not isinstance(observed, desired):
                return c_exc.UseCaseError(
                    f"argument `{observed}` should be a instance of `{desired}`.",
                    exp=True,
                    logger=LOGGER,
                )()

        # ? --------------------------------------------------------------------
        # ? Calculate joint probability units
        # ? --------------------------------------------------------------------

        group: LabeledPriors
        clade_adherence_stats: dict[PriorGroup, float] = {}

        if isinstance(clade_priors, IngroupCladePriors):
            input_labeled_priors: IngroupLabeledPriors = next(
                i for i in clade_priors.priors if i.group == PriorGroup.INGROUP
            )

            sister_labeled_priors: SisterGroupLabeledPriors = next(
                i for i in clade_priors.priors if i.group == PriorGroup.SISTER
            )

            if calculate_bootstrap is True:
                noise_group_priors = {
                    **input_labeled_priors.priors,
                    **sister_labeled_priors.priors,
                }

                noise_labels = set(
                    sorted(
                        [
                            *list(input_labeled_priors.labels),
                            *list(sister_labeled_priors.labels),
                        ]
                    )
                )

            for group in [
                input_labeled_priors,
                sister_labeled_priors,
            ]:
                if (
                    adherence_either := calculate_clade_adherence_test(
                        labeled_priors=group,
                        query_kmers=query_kmers,
                        kmer_indices=kmer_indices,
                        total_length=len(group.labels),
                        adherence_strategy=adherence_strategy,
                    )
                ).is_left:
                    return adherence_either

                adherence = adherence_either.value

                if calculate_bootstrap is True:
                    if (
                        calculation_either := __calculate_bootstrap(
                            query_kmers=query_kmers,
                            kmer_indices=kmer_indices,
                            bootstrap_replicates=bootstrap_replicates,
                            sample_size=len(group.priors),
                            noise_group_priors=copy(noise_group_priors),
                            noise_group_labels=copy(noise_labels),
                            adherence_strategy=adherence_strategy,
                        )
                    ).is_left:
                        return calculation_either

                    (
                        bootstrap_joint_probabilities,
                        bootstrap_matching_kmers,
                    ) = calculation_either.value

                    if (
                        adherence_strategy
                        == AdherenceTestStrategy.JOINT_PROBABILITY
                    ):
                        adherence.adherence_p_value = (
                            sum(
                                1
                                for i in bootstrap_joint_probabilities
                                if i <= adherence.joint_probability
                            )
                            / bootstrap_replicates
                        )

                    else:
                        adherence.adherence_p_value = (
                            sum(
                                1
                                for i in bootstrap_matching_kmers
                                if i >= adherence.match_kmers
                            )
                            / bootstrap_replicates
                        )

                LOGGER.debug(f"\t\t{group.group}: {adherence}")

                clade_adherence_stats.update({group.group: adherence})

        if isinstance(clade_priors, OutgroupCladePriors):
            if (
                adherence_either := calculate_clade_adherence_test(
                    labeled_priors=clade_priors.priors,
                    query_kmers=query_kmers,
                    kmer_indices=kmer_indices,
                    total_length=len(clade_priors.priors.labels),
                    adherence_strategy=adherence_strategy,
                )
            ).is_left:
                return adherence_either

            clade_adherence_stats.update(
                {clade_priors.priors.group: adherence_either.value}
            )

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(clade_adherence_stats)

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()


def __calculate_bootstrap(
    query_kmers: set[str],
    kmer_indices: KmersInverseIndices,
    bootstrap_replicates: int,
    sample_size: int,
    noise_group_priors: dict[str, float],
    noise_group_labels: set[str],
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
                labels=noise_group_labels,
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
