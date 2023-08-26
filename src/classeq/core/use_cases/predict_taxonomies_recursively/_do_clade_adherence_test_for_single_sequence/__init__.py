from random import seed

import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.priors import (
    CladePriors,
    IngroupCladePriors,
    OutgroupCladePriors,
    PriorGroup,
)
from classeq.settings import LOGGER

from ._dtos import AdherenceResult, AdherenceStatus


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

    def measure_intersection(
        clade_priors: CladePriors,
    ) -> AdherenceResult:
        """Measure the intersection of the query kmers and the priors.

        Attempts to avoid change the `query_kmers` variable which is inherited
        from the parent function context.

        Args:
            clade_priors (LabeledPriors): The labeled priors.

        Returns:
            AdherenceResult: The adherence result.

        """

        prior_keys = clade_priors.kmers
        match_kmers = query_kmers.intersection(set(prior_keys)).__len__()
        status = AdherenceStatus.UNDEFINED

        if match_kmers == 0:
            status = AdherenceStatus.NOT_ENOUGH_PRIORS
        else:
            status = AdherenceStatus.SUCCESS

        return AdherenceResult(
            match_kmers=match_kmers,
            query_kmers_size=len(query_kmers),
            subject_kmers_size=len(prior_keys),
            status=status,
        )

    try:
        seed(987_654_321)

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

        labeled_prior: CladePriors
        clade_adherence_stats: dict[PriorGroup, float] = {}

        # ? --------------------------------------------------------------------
        # ? Calculate adherence test for Outgroup clades
        # ? --------------------------------------------------------------------

        if isinstance(clade_priors, OutgroupCladePriors):
            adherence_result = measure_intersection(
                clade_priors=clade_priors.clade_priors
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

        if isinstance(clade_priors, IngroupCladePriors):
            for labeled_prior in clade_priors.clade_priors:
                adherence_result = measure_intersection(
                    clade_priors=labeled_prior
                )

                LOGGER.debug(f"\t\t{labeled_prior.group}: {adherence_result}")

                clade_adherence_stats.update(
                    {labeled_prior.group: adherence_result}
                )

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(clade_adherence_stats)

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
