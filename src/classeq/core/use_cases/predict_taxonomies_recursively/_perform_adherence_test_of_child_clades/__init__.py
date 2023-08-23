from typing import Any

import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.domain.dtos.clade import ClasseqClade
from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.priors import PriorGroup, TreePriors
from classeq.settings import LOGGER

from .._calculate_clade_adherence_with_bootstrap._dtos import (
    AdherenceResult,
    AdherenceTestStrategy,
)
from .._do_clade_adherence_test_for_single_sequence import (
    do_clade_adherence_test_for_single_sequence,
)
from ._dtos import CladeAdherenceResult, CladeAdherenceResultStatus


def perform_adherence_test_of_child_clades(
    query_kmers: set[str],
    clades: list[ClasseqClade],
    tree_priors: TreePriors,
    kmer_indices: KmersInverseIndices,
    **kwargs: Any,
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
        sister_joint_probabilities: list[AdherenceResult] = []

        adherence_strategy: AdherenceTestStrategy = kwargs.get(
            "adherence_strategy", AdherenceTestStrategy(None)
        )

        for clade in clades:
            LOGGER.debug("")
            LOGGER.debug(f"\tProcessing clade: {clade}")

            try:
                clade_priors = next(
                    i for i in tree_priors.ingroups if i.parent == clade.id
                )
            except StopIteration:
                LOGGER.warning(f"Ignore child: {clade}")
                continue

            if (
                adherence_test_either := do_clade_adherence_test_for_single_sequence(
                    query_kmers=query_kmers,
                    clade_priors=clade_priors,
                    kmer_indices=kmer_indices,
                    **kwargs,
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

            if adherence_strategy == AdherenceTestStrategy.JOINT_PROBABILITY:
                if ingroup.joint_probability < sister.joint_probability:
                    contrasting_clades.add(
                        CladeAdherenceResult(
                            clade=clade,
                            ingroup_adherence_test=ingroup,
                            sister_adherence_test=sister,
                        )
                    )
            else:
                if ingroup.match_kmers > sister.match_kmers:
                    contrasting_clades.add(
                        CladeAdherenceResult(
                            clade=clade,
                            ingroup_adherence_test=ingroup,
                            sister_adherence_test=sister,
                        )
                    )

            sister_joint_probabilities.append(sister)

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        if len(contrasting_clades) == 0:
            return right(
                (
                    (
                        min(
                            sister_joint_probabilities,
                            key=lambda i: i.joint_probability,
                        )
                        if len(sister_joint_probabilities) > 0
                        else None
                    ),
                    None,
                    CladeAdherenceResultStatus.MAX_RESOLUTION_REACHED,
                )
            )

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

        if len(contrasting_clades) > 1:
            return right(
                (
                    sorted(
                        contrasting_clades,
                        key=lambda i: i.ingroup_adherence_test.joint_probability,
                    )[0].ingroup_adherence_test,
                    None,
                    CladeAdherenceResultStatus.INCONCLUSIVE,
                )
            )

        return c_exc.UseCaseError(
            "Unable to determine the clade adherence test result.",
            logger=LOGGER,
        )()

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
