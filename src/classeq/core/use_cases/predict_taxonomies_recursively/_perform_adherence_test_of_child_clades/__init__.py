from copy import deepcopy

import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.domain.dtos.clade import CladeWrapper
from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.priors import PriorGroup, TreePriors
from classeq.settings import LOGGER

from .._do_clade_adherence_test_for_single_sequence import (
    do_clade_adherence_test_for_single_sequence,
)
from ._dtos import CladeAdherenceResult, CladeAdherenceResultStatus


def perform_adherence_test_of_child_clades(
    query_kmers: list[str],
    clades: list[CladeWrapper],
    tree_priors: TreePriors,
    kmer_indices: KmersInverseIndices,
) -> Either[
    c_exc.MappedErrors,
    tuple[float, CladeAdherenceResult | None, CladeAdherenceResultStatus],
]:
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
        ingroup_joint_probabilities: list[float] = []
        sister_joint_probabilities: list[float] = []

        for clade in clades:
            try:
                clade_priors = next(
                    i for i in tree_priors.ingroups if i.parent == clade.id
                )
            except StopIteration:
                LOGGER.debug(f"Ignore child: {clade}")
                raise NotImplementedError()

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

            if ingroup < sister:
                contrasting_clades.add(
                    CladeAdherenceResult(
                        clade=deepcopy(clade),
                        ingroup_joint_probability=ingroup,
                        sister_joint_probability=sister,
                    )
                )

            ingroup_joint_probabilities.append(ingroup)
            sister_joint_probabilities.append(sister)

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        if len(contrasting_clades) == 0:
            return right(
                (
                    min(sister_joint_probabilities),
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
                        clade_binding.ingroup_joint_probability,
                        clade_binding,
                        CladeAdherenceResultStatus.NEXT_ITERATION,
                    )
                )

            return right(
                (
                    clade_binding.ingroup_joint_probability,
                    clade_binding,
                    CladeAdherenceResultStatus.CONCLUSIVE_INGROUP,
                )
            )

        if len(contrasting_clades) > 1:
            return right(
                (
                    sorted(
                        contrasting_clades,
                        key=lambda i: i.ingroup_joint_probability,
                    )[0],
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
