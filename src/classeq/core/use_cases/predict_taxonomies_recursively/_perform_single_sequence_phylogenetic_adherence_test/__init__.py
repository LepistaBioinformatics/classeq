from copy import copy
from typing import Any, Literal

import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.domain.dtos.clade import ClasseqClade
from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.priors import (
    OutgroupPriors,
    PriorGroup,
    TreePriors,
)
from classeq.core.domain.dtos.strand import StrandEnum
from classeq.settings import DEFAULT_MATCHES_COVERAGE, LOGGER

from .._do_clade_adherence_test_for_single_sequence import (
    do_clade_adherence_test_for_single_sequence,
)
from .._do_clade_adherence_test_for_single_sequence._dtos import (
    AdherenceResult,
    AdherenceStatus,
)
from .._perform_adherence_test_of_child_clades import (
    CladeAdherenceResult,
    CladeAdherenceResultStatus,
    perform_adherence_test_of_child_clades,
)
from ._dtos import PredictionStep, PredictionResult


def perform_single_sequence_phylogenetic_adherence_test(
    target_sequence: str,
    kmers_indices: KmersInverseIndices,
    tree_priors: TreePriors,
    kmer_size: int,
    tree: ClasseqClade,
    strand: StrandEnum,
    outgroup_priors: OutgroupPriors | Literal[False] = False,
    max_iterations: int = 1000,
    matches_coverage: float = DEFAULT_MATCHES_COVERAGE,
    **kwargs: Any,
) -> Either[c_exc.MappedErrors, PredictionResult]:
    """Perform phylogenetic adherence test.

    Description:
        This function performs phylogenetic adherence test. The test is
        performed by recursively traversing the tree and performing
        clade-adherence test for each clade. The test is performed until
        the test is conclusive or inconclusive.

    Args:
        target_sequence (str): The target sequence.
        kmers_indices (KmersInverseIndices): The kmer indices.
        tree_priors (TreePriors): The tree priors.
        kmer_size (int): The kmer size.
        tree (ClasseqClade): The tree.
        strand (StrandEnum): The strand.
        max_iterations (int, optional): The maximum number of iterations.
            Defaults to 1000.
        matches_coverage (float, optional): The minimum number of matches
            coverage. Defaults to 0.5.

    Returns:
        Either[c_exc.MappedErrors, bool]: The result of the test.

    Raises:
        c_exc.UseCaseError: If the retrieved tree using `get_hierarchical_tree`
            is not a rooted tree.

    """

    try:
        # ? --------------------------------------------------------------------
        # ? Validate args
        # ? --------------------------------------------------------------------

        if not isinstance(target_sequence, str):
            return c_exc.UseCaseError(
                "Unexpected error. The `target_sequence` argument should be a "
                + f"string. Received {type(target_sequence)}.",
                logger=LOGGER,
            )()

        if not isinstance(tree_priors, TreePriors):
            return c_exc.UseCaseError(
                "Unexpected error. The `tree_priors` argument should be a "
                + f"TreePriors. Received {type(tree_priors)}.",
                logger=LOGGER,
            )()

        if max_iterations > 9999:
            return c_exc.UseCaseError(
                "Unexpected error. The `max_iterations` argument should be "
                + f"less than 9999. Received {max_iterations}.",
                logger=LOGGER,
            )()

        if tree.is_root() is False:
            return c_exc.UseCaseError(
                "Unexpected error. Retrieved tree using "
                + "`get_hierarchical_tree` is not a rooted tree.",
                logger=LOGGER,
            )()

        if not isinstance(kmers_indices, KmersInverseIndices):
            return c_exc.UseCaseError(
                "Unexpected error. The `kmers_indices` argument should be a "
                + f"KmersInverseIndices. Received {type(kmers_indices)}.",
                logger=LOGGER,
            )()

        if not isinstance(kmer_size, int):
            return c_exc.UseCaseError(
                "Unexpected error. The `kmer_size` argument should be a "
                + f"int. Received {type(kmer_size)}.",
                logger=LOGGER,
            )()

        if not isinstance(tree, ClasseqClade):
            return c_exc.UseCaseError(
                "Unexpected error. The `tree` argument should be a "
                + f"ClasseqClade. Received {type(tree)}.",
                logger=LOGGER,
            )()

        if not isinstance(strand, StrandEnum):
            return c_exc.UseCaseError(
                "Unexpected error. The `strand` argument should be a "
                + f"StrandEnum. Received {type(strand)}.",
                logger=LOGGER,
            )()

        if matches_coverage < 0.00001 or matches_coverage > 1:
            return c_exc.UseCaseError(
                "Unexpected error. The `matches_coverage` argument should be "
                + f"between 0.00001 and 1. Received {matches_coverage}.",
                logger=LOGGER,
            )()

        if outgroup_priors:
            if not isinstance(outgroup_priors, OutgroupPriors):
                return c_exc.UseCaseError(
                    "Unexpected error. The `outgroup_priors` argument should be a "
                    + f"OutgroupPriors. Received {type(outgroup_priors)}.",
                    logger=LOGGER,
                )()

        # ? --------------------------------------------------------------------
        # ? Get sequence kmers
        # ? --------------------------------------------------------------------

        query_kmers: set[str] = {
            kmer
            for kmer in KmersInverseIndices.generate_kmers(
                dna_sequence=KmersInverseIndices.sanitize_sequence(
                    target_sequence,
                    as_seq=True,
                ),
                k_size=kmer_size,
                strand=strand,
            )
        }

        # ? --------------------------------------------------------------------
        # ? Perform adherence test for outgroup
        # ? --------------------------------------------------------------------

        if outgroup_priors:
            if (
                binding_either := do_clade_adherence_test_for_single_sequence(
                    query_kmers=query_kmers,
                    clade_priors=outgroup_priors,
                    kmer_indices=kmers_indices,
                )
            ).is_left:
                return binding_either

            if (
                outgroup_adherence_test := binding_either.value.pop(
                    PriorGroup.OUTGROUP
                )
            ) is None:
                return c_exc.UseCaseError(
                    "Unexpected error on try to calculate "
                    + "adherence test for ingroup.",
                    logger=LOGGER,
                )()

            LOGGER.debug(f"Output Adherence Result: {outgroup_adherence_test}")

        # ? --------------------------------------------------------------------
        # ? Perform adherence test for the outgroup-sister pairs
        # ? --------------------------------------------------------------------

        if (ingroup_clades_either := tree.get_ingroup_clades()).is_left:
            return ingroup_clades_either

        ingroup_clades: list[ClasseqClade] = list(
            ingroup for ingroup in ingroup_clades_either.value
        )

        # ? --------------------------------------------------------------------
        # ? Perform adherence test for the ingroups
        # ? --------------------------------------------------------------------

        response_children: CladeAdherenceResult | None
        adherence_result: AdherenceResult
        status = CladeAdherenceResultStatus.NEXT_ITERATION
        local_max_iterations = copy(max_iterations)
        response_clades: list[ClasseqClade] = ingroup_clades
        current_iteration: int = 0

        clade_path: list[PredictionStep] = list()

        while (
            response_clades
            and status == CladeAdherenceResultStatus.NEXT_ITERATION
        ):
            current_iteration += 1

            if len(response_clades) == 1:
                clade = response_clades.pop(0)
                if (children := clade.children) is None:
                    continue
            else:
                children = response_clades
                response_clades = list()

            children = [i for i in children if i.is_internal()]

            if len(children) == 0:
                break

            if (
                response_either := perform_adherence_test_of_child_clades(
                    query_kmers=query_kmers,
                    clades=children,
                    tree_priors=tree_priors,
                    kmer_indices=kmers_indices,
                    **kwargs,
                )
            ).is_left:
                return response_either

            (
                adherence_result,
                response_children,
                status,
            ) = response_either.value

            LOGGER.debug("")
            LOGGER.debug(f"\tChildren Adherence Result: {response_children}")
            LOGGER.debug("")
            LOGGER.debug(f"\tStatus: {status}")

            if response_children is not None:
                clade_path.append(
                    PredictionStep(
                        depth=current_iteration,
                        result=response_children,
                    )
                )

            # ------------------------------------------------------------------
            # Break the search if the adherence test is inconclusive due to
            # absence of priors to perform comparisons.
            # ------------------------------------------------------------------
            if status in [
                CladeAdherenceResultStatus.MAX_RESOLUTION_REACHED,
                CladeAdherenceResultStatus.CONCLUSIVE_INGROUP,
            ]:
                break

            # ------------------------------------------------------------------
            # Break the search if the adherence test is inconclusive due to
            # absence of priors to perform comparisons.
            # ------------------------------------------------------------------
            if adherence_result.status == AdherenceStatus.NOT_ENOUGH_PRIORS:
                LOGGER.warning(
                    f"Unable to classify clade `{'clade'}`. Not enough "
                    + "priors evidences to run comparisons."
                )

                status = CladeAdherenceResultStatus(None)
                break

            # ------------------------------------------------------------------
            # Break the search if the first iteration is inconclusive for
            # ingroup. The first iteration represents the test for the sister
            # clade of the outgroup, in other words, the tree entrypoint.
            # ------------------------------------------------------------------

            if outgroup_priors:
                if (
                    current_iteration == 1
                    and outgroup_adherence_test.match_kmers
                    > adherence_result.match_kmers
                ):
                    LOGGER.debug(
                        "The processed sequence does not differs from outgroup"
                    )

                    clade_path[0].result.clade = ClasseqClade.new_anemic_clade()
                    status = CladeAdherenceResultStatus.CONCLUSIVE_OUTGROUP
                    break

            # ------------------------------------------------------------------
            # Case the adherence test is conclusive for ingroup and the ingroup
            # clade has children clades, the adherence test should return a
            # `CladeAdherenceResultStatus.NEXT_ITERATION`, and the children of
            # the ingroup clade should be added to the list of clades to be
            # testes.
            # ------------------------------------------------------------------
            if (
                status == CladeAdherenceResultStatus.NEXT_ITERATION
                and response_children is not None
            ):
                # --------------------------------------------------------------
                # Break the search if the number of match kmers not exceeds the
                # at last 50% of the total number of query kmers.
                # --------------------------------------------------------------

                if adherence_result.match_kmers < (
                    query_kmers.__len__() * matches_coverage
                ):
                    LOGGER.debug(
                        "Low number of matches. The processed sequence "
                        + "probably does not belongs to the current clade."
                    )

                    clade_path[
                        current_iteration - 1
                    ].result.clade = ClasseqClade.new_anemic_clade()

                    status = CladeAdherenceResultStatus.INCONCLUSIVE
                    break

                response_clades.append(response_children.clade)

            # ------------------------------------------------------------------
            # Otherwise, the adherence test is conclusive the search loop should
            # be broken.
            # ------------------------------------------------------------------
            else:
                break

            # ------------------------------------------------------------------
            # This is a safety measure to avoid infinite loops. Don't remove it.
            # ------------------------------------------------------------------
            local_max_iterations -= 1
            if local_max_iterations == 0:
                LOGGER.critical("Max iterations reached.")
                break

        if status == CladeAdherenceResultStatus.NEXT_ITERATION:
            status = CladeAdherenceResultStatus.CONCLUSIVE_INGROUP
            LOGGER.debug("Max resolution reached")

        if clade_path.__len__() == 0:
            status = CladeAdherenceResultStatus(None)

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(
            PredictionResult(
                status=status,
                path=clade_path,
            )
        )

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
