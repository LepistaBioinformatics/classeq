import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.ordered_tuple import OrderedTuple
from classeq.core.domain.dtos.priors import (
    CladePriors,
    IngroupPriors,
    OutgroupPriors,
    PriorGroup,
)
from classeq.settings import LOGGER

from ._dtos import AdherenceResult, AdherenceStatus


def do_clade_adherence_test_for_single_sequence(
    query_kmers: set[str],
    clade_priors: IngroupPriors | OutgroupPriors,
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

    def measure_intersection(kmers: OrderedTuple) -> AdherenceResult:
        """Measure the intersection of the query kmers and the priors.

        Attempts to avoid change the `query_kmers` variable which is inherited
        from the parent function context.

        Args:
            clade_priors (LabeledPriors): The labeled priors.

        Returns:
            AdherenceResult: The adherence result.

        """

        match_kmers = query_kmers.intersection(set(kmers)).__len__()

        status = AdherenceStatus.UNDEFINED
        if match_kmers == 0:
            status = AdherenceStatus.NOT_ENOUGH_PRIORS
        else:
            status = AdherenceStatus.SUCCESS

        return AdherenceResult(
            match_kmers=match_kmers,
            query_kmers_size=query_kmers.__len__(),
            subject_kmers_size=kmers.__len__(),
            status=status,
        )

    try:
        # ? --------------------------------------------------------------------
        # ? Validate entries
        # ? --------------------------------------------------------------------

        if not any(
            [
                isinstance(clade_priors, IngroupPriors),
                isinstance(clade_priors, OutgroupPriors),
            ]
        ):
            return c_exc.UseCaseError(
                f"argument `{clade_priors}` should be any of "
                + f"`{IngroupPriors}` or `{OutgroupPriors}`.",
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

        target_clade_priors: CladePriors
        clade_adherence_stats: dict[PriorGroup, float] = {}

        # ? --------------------------------------------------------------------
        # ? Calculate adherence test for Outgroup clades
        # ? --------------------------------------------------------------------

        if isinstance(clade_priors, OutgroupPriors):
            adherence_result = measure_intersection(
                kmers=clade_priors.clade_priors.kmers
            )

            LOGGER.debug(
                f"\t\t{clade_priors.clade_priors.group}: {adherence_result}"
            )

            clade_adherence_stats.update(
                {clade_priors.clade_priors.group: adherence_result}
            )

        # ? --------------------------------------------------------------------
        # ? Calculate adherence test for Ingroup clades
        # ? --------------------------------------------------------------------

        if isinstance(clade_priors, IngroupPriors):
            for target_clade_priors in clade_priors.clade_priors:
                adherence_result = measure_intersection(
                    kmers=target_clade_priors.kmers
                )

                LOGGER.debug(
                    f"\t\t{target_clade_priors.group}: {adherence_result}"
                )

                clade_adherence_stats.update(
                    {target_clade_priors.group: adherence_result}
                )

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(clade_adherence_stats)

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
