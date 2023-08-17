import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.priors import (
    IngroupCladePriors,
    OutgroupCladePriors,
    PriorGroup,
)
from classeq.settings import LOGGER

from ._calculate_clade_adherence_with_bootstrap import (
    AdherenceResult,
    calculate_clade_adherence_with_bootstrap,
)


def do_clade_adherence_test_for_single_sequence(
    query_kmers: set[str],
    clade_priors: IngroupCladePriors | OutgroupCladePriors,
    kmer_indices: KmersInverseIndices,
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

        clade_adherence_stats: dict[PriorGroup, float] = {}

        if isinstance(clade_priors, IngroupCladePriors):
            for group in sorted(
                clade_priors.priors,
                key=lambda x: x.group.value,
            ):
                if (
                    adherence_either := calculate_clade_adherence_with_bootstrap(
                        labeled_priors=group,
                        query_kmers=query_kmers,
                        kmer_indices=kmer_indices,
                        total_length=len(group.labels),
                    )
                ).is_left:
                    return adherence_either

                LOGGER.debug(f"\t\t{group.group}: {adherence_either.value}")

                clade_adherence_stats.update(
                    {group.group: adherence_either.value}
                )

        if isinstance(clade_priors, OutgroupCladePriors):
            if (
                adherence_either := calculate_clade_adherence_with_bootstrap(
                    labeled_priors=clade_priors.priors,
                    query_kmers=query_kmers,
                    kmer_indices=kmer_indices,
                    total_length=len(clade_priors.priors.labels),
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
