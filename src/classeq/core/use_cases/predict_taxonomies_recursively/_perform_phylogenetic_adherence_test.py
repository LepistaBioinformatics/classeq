from copy import copy, deepcopy

import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.domain.dtos.clade import CladeWrapper
from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.priors import PriorGroup, TreePriors
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.settings import DEFAULT_KMER_SIZE, LOGGER

from ._do_clade_adherence_test_for_single_sequence import (
    do_clade_adherence_test_for_single_sequence,
)
from ._perform_adherence_test_of_child_clades import (
    CladeAdherenceResult,
    CladeAdherenceResultStatus,
    perform_adherence_test_of_child_clades,
)


def perform_phylogenetic_adherence_test(
    target_sequence: str,
    reference_set: ReferenceSet,
    tree_priors: TreePriors,
    max_iterations: int = 1000,
) -> Either[c_exc.MappedErrors, CladeWrapper | None]:
    """Perform phylogenetic adherence test.

    Description:
        This function performs phylogenetic adherence test. The test is
        performed by recursively traversing the tree and performing
        clade-adherence test for each clade. The test is performed until
        the test is conclusive or inconclusive.

    Args:
        target_sequence (str): The sequence to be tested.
        reference_set (ReferenceSet): The reference set.
        tree_priors (TreePriors): The tree priors.

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

        if not isinstance(reference_set, ReferenceSet):
            return c_exc.UseCaseError(
                "Unexpected error. The `reference_set` argument should be a "
                + f"ReferenceSet. Received {type(reference_set)}.",
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

        # ? --------------------------------------------------------------------
        # ? Get sequence kmers
        # ? --------------------------------------------------------------------

        target_sequence_kmers: list[str] = [
            kmer
            for kmer in KmersInverseIndices.generate_kmers(
                dna_sequence=target_sequence.upper(),
                k_size=DEFAULT_KMER_SIZE,
            )
        ]

        # ? --------------------------------------------------------------------
        # ? Collect hierarchical tree
        # ? --------------------------------------------------------------------

        if (tree_either := reference_set.get_hierarchical_tree()).is_left:
            return tree_either

        tree: CladeWrapper = tree_either.value

        if tree.is_root() is False:
            return c_exc.UseCaseError(
                "Unexpected error. Retrieved tree using "
                + "`get_hierarchical_tree` is not a rooted tree.",
                logger=LOGGER,
            )()

        # ? --------------------------------------------------------------------
        # ? Perform adherence test for outgroup
        # ? --------------------------------------------------------------------

        if (outgroup_clades_either := tree.get_outgroup_clade()).is_left:
            return outgroup_clades_either

        outgroup_clades: CladeWrapper = outgroup_clades_either.value

        if (
            binding_either := do_clade_adherence_test_for_single_sequence(
                query_kmers=target_sequence_kmers,
                clade_priors=tree_priors.outgroup,
                kmer_indices=reference_set.msa.kmers_indices,
            )
        ).is_left:
            return binding_either

        if (
            outgroup_joint_probability := binding_either.value.pop(
                PriorGroup.OUTGROUP
            )
        ) is None:
            return c_exc.UseCaseError(
                "Unexpected error on try to calculate "
                + "adherence test for ingroup.",
                logger=LOGGER,
            )()

        # ? --------------------------------------------------------------------
        # ? Perform adherence test for the outgroup-sister pairs
        # ? --------------------------------------------------------------------

        if (ingroup_clades_either := tree.get_ingroup_clades()).is_left:
            return ingroup_clades_either

        ingroup_clades: list[CladeWrapper] = list(
            ingroup
            for ingroup in ingroup_clades_either.value
            if ingroup.parent == tree_priors.outgroup.parent
        )

        # ? --------------------------------------------------------------------
        # ? Perform adherence test for the ingroups
        # ? --------------------------------------------------------------------

        response_children: CladeAdherenceResult | None
        joint_probability: float
        status = CladeAdherenceResultStatus.NEXT_ITERATION
        local_max_iterations = copy(max_iterations)
        final_response: CladeWrapper | None = None
        response_clades: list[CladeWrapper] = deepcopy(ingroup_clades)
        current_iteration: int = 0

        while (
            response_clades
            and status == CladeAdherenceResultStatus.NEXT_ITERATION
        ):
            clade = response_clades.pop(0)
            current_iteration += 1

            if (children := clade.children) is None:
                continue

            children = [i for i in children if i.is_internal()]

            if len(children) == 0:
                final_response = clade
                break

            if (
                response_either := perform_adherence_test_of_child_clades(
                    query_kmers=target_sequence_kmers,
                    clades=children,
                    tree_priors=tree_priors,
                    kmer_indices=reference_set.msa.kmers_indices,
                )
            ).is_left:
                return response_either

            (
                joint_probability,
                response_children,
                status,
            ) = response_either.value

            # ------------------------------------------------------------------
            # Break the search if the first iteration is inconclusive for
            # ingroup. The first iteration represents the test for the sister
            # clade of the outgroup, in other words, the tree entrypoint.
            # ------------------------------------------------------------------
            if (
                current_iteration == 1
                and outgroup_joint_probability < joint_probability
            ):
                LOGGER.debug("The processed sequence is an outgroup.")
                final_response = outgroup_clades
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
                response_clades.append(response_children.clade)

            # ------------------------------------------------------------------
            # Otherwise, the adherence test is conclusive the search loop should
            # be broken.
            # ------------------------------------------------------------------
            else:
                final_response = clade
                break

            # ------------------------------------------------------------------
            # This is a safety measure to avoid infinite loops. Don't remove it.
            # ------------------------------------------------------------------
            local_max_iterations -= 1
            if local_max_iterations == 0:
                LOGGER.critical("Max iterations reached.")
                break

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(final_response)

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
