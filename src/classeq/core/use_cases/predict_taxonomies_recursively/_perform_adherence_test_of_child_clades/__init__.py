from typing import Any

import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.domain.dtos.clade import ClasseqClade
from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.priors import PriorGroup, TreePriors
from classeq.settings import (
    LOGGER,
    MINIMUM_INGROUP_QUERY_KMERS_MATCH,
    MINIMUM_INGROUP_SISTER_MATCH_KMERS_DIFFERENCE,
)

from .._do_clade_adherence_test_for_single_sequence import (
    do_clade_adherence_test_for_single_sequence,
    AdherenceStatus,
)
from .._do_clade_adherence_test_for_single_sequence._dtos import AdherenceResult
from ._dtos import CladeAdherenceResult, CladeAdherenceResultStatus


def perform_adherence_test_of_child_clades(
    query_kmers: set[str],
    clades: list[ClasseqClade],
    tree_priors: TreePriors,
    kmer_indices: KmersInverseIndices,
    minimum_ingroup_query_kmers_match: int = MINIMUM_INGROUP_QUERY_KMERS_MATCH,
    ingroup_sister_match_kmers_difference: int = MINIMUM_INGROUP_SISTER_MATCH_KMERS_DIFFERENCE,
    **_: Any,
) -> Either[
    c_exc.MappedErrors,
    tuple[
        AdherenceResult, CladeAdherenceResult | None, CladeAdherenceResultStatus
    ],
]:
    """Perform adherence test of child clades.

    Description:
        This function performs adherence test of child clades. The adherence
        test is performed to all child clades and the results are stored to
        further comparisons. Given the possibility to more than one clade would
        selected all clades should be tested and the output results further
        evaluated.

    Args:
        query_kmers (set[str]): The query kmers.
        clades (list[CladeWrapper]): The clades.
        tree_priors (TreePriors): The tree priors.
        kmer_indices (KmersInverseIndices): The kmer indices.
        minimum_query_kmers_match (int, optional): The minimum query kmers
            match. Defaults to 10.
        ingroup_sister_match_kmers_difference (int, optional): The ingroup
            sister match kmers difference. Defaults to 10.

    Returns:
        Either[
            c_exc.MappedErrors,
            tuple[
                AdherenceResult,
                CladeAdherenceResult | None,
                CladeAdherenceResultStatus,
            ],
        ]: The adherence test result.

    Raises:
        c_exc.UseCaseError: If any of the arguments is not a list.
        c_exc.UseCaseError: If any of the sub-items of the arguments is
            not a list.

    """

    LOGGER.debug("")
    LOGGER.debug("perform_adherence_test_of_child_clades")

    try:
        if len(clades) == 0:
            return c_exc.UseCaseError(
                "Clades should be an array with at last one element.",
                logger=LOGGER,
            )()

        # ? --------------------------------------------------------------------
        # ? Run adherence test for all clades
        #
        # Adherence test should be executed to all input clades and results
        # stored to further comparisons. Given the possibility to more than one
        # clade would selected all clades should be tested and the output
        # results further evaluated.
        #
        # ? --------------------------------------------------------------------

        contrasting_clades: set[CladeAdherenceResult] = set()
        sister_adherence_tests: list[AdherenceResult] = []

        for clade in clades:
            LOGGER.debug("")
            LOGGER.debug(f"\tProcessing clade: {clade}")

            try:
                clade_priors = next(
                    i for i in tree_priors.ingroups if i.parent == clade.id
                )
            except StopIteration:
                LOGGER.debug(f"Ignore child with by small clade size: {clade}")
                continue

            if (
                adherence_test_either := do_clade_adherence_test_for_single_sequence(
                    query_kmers=query_kmers,
                    clade_priors=clade_priors,
                    kmer_indices=kmer_indices,
                )
            ).is_left:
                return adherence_test_either

            adherence_test = adherence_test_either.value

            if (ingroup := adherence_test.pop(PriorGroup.INGROUP)) is None:
                return c_exc.UseCaseError(
                    "Unexpected error on try to calculate "
                    + "adherence test for ingroup.",
                    logger=LOGGER,
                )()

            if (sister := adherence_test.pop(PriorGroup.SISTER)) is None:
                return c_exc.UseCaseError(
                    "Unexpected error on try to calculate "
                    + "adherence test for sister.",
                    logger=LOGGER,
                )()

            # ? ----------------------------------------------------------------
            # ? Evaluate adherence test results
            #
            # Adherence test results should be evaluated to determine if the
            # clade is a good candidate to be selected as the best clade. These
            # step is performed by comparing the adherence test results of the
            # ingroup and sister group. The following conditions should be
            # satisfied to consider the clade as a good candidate:
            #
            # ? ----------------------------------------------------------------

            if all(
                [
                    #
                    # 1. Ingroup number of match kmers should be greater than
                    # sister group matches.
                    #
                    ingroup.match_kmers > sister.match_kmers,
                    #
                    # 2. Ingroup and sister group kmers match evaluation should
                    # occurred with success.
                    #
                    (
                        ingroup.status == AdherenceStatus.SUCCESS
                        and sister.status == AdherenceStatus.SUCCESS
                    ),
                    #
                    # 3. Ingroup number of match kmers should be greater than
                    # minimum query kmers match. Otherwise the difference
                    # between ingroup and sister should be considered random.
                    #
                    ingroup.match_kmers >= minimum_ingroup_query_kmers_match,
                    #
                    # 4. Ingroup and sister group match kmers difference should
                    # be greater than minimum query kmers match. Otherwise the
                    # difference between ingroup and sister should be considered
                    # random.
                    #
                    (
                        (ingroup.match_kmers - sister.match_kmers)
                        >= ingroup_sister_match_kmers_difference
                    ),
                ]
            ):
                contrasting_clades.add(
                    CladeAdherenceResult(
                        clade=clade,
                        ingroup_adherence_test=ingroup,
                        sister_adherence_test=sister,
                    )
                )

            sister_adherence_tests.append(sister)

        # ? --------------------------------------------------------------------
        # ? Evaluate selected best clades
        #
        # Clades selected (or not) at the previous step should be evaluated to
        # build the use-case output response.
        #
        # ? --------------------------------------------------------------------
        #
        # If no clade was selected, the best results of the sister clade will be
        # used to build the final response. In this situation the response
        # status should be `MAX_RESOLUTION_REACHED`, indicating that the maximum
        # resolution inside the tree clades was reached.
        #
        if len(contrasting_clades) == 0:
            return right(
                (
                    (
                        max(
                            sister_adherence_tests,
                            key=lambda i: i.match_kmers,
                        )
                        if len(sister_adherence_tests) > 0
                        else None
                    ),
                    None,
                    CladeAdherenceResultStatus.MAX_RESOLUTION_REACHED,
                )
            )
        #
        # If only one clade was selected, final results is further evaluated to
        # determine if the clade is conclusive (CONCLUSIVE_INGROUP) or the tree
        # search should continue (NEXT_ITERATION).
        #
        if len(contrasting_clades) == 1:
            clade_binding: CladeAdherenceResult

            if (clade_binding := contrasting_clades.pop()) is None:
                return c_exc.UseCaseError(
                    "Unexpected error on try to get clade from set.",
                    logger=LOGGER,
                )()

            if clade_binding.clade.children is not None:
                return right(
                    (
                        clade_binding.ingroup_adherence_test,
                        clade_binding,
                        CladeAdherenceResultStatus.NEXT_ITERATION,
                    )
                )

            return right(
                (
                    clade_binding.ingroup_adherence_test,
                    clade_binding,
                    CladeAdherenceResultStatus.CONCLUSIVE_INGROUP,
                )
            )
        #
        # If more than one clade was selected, the clade with the highest
        # number of match kmers will be selected. In this situation the response
        # status should be `INCONCLUSIVE` (more than one clade was selected).
        #
        if len(contrasting_clades) > 1:
            return right(
                (
                    sorted(
                        contrasting_clades,
                        key=lambda i: i.ingroup_adherence_test.match_kmers,
                        reverse=True,
                    )[0].ingroup_adherence_test,
                    None,
                    CladeAdherenceResultStatus(None),
                )
            )
        #
        # If none of the previous conditions was satisfied, an unexpected error
        # should be raised.
        #
        return c_exc.UseCaseError(
            "Unable to determine the clade adherence test result.",
            logger=LOGGER,
        )()

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
