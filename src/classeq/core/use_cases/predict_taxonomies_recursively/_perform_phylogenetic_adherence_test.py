import clean_base.exceptions as c_exc
from clean_base.either import Either, left, right

from classeq.core.domain.dtos.clade import CladeWrapper
from classeq.core.domain.dtos.priors import TreePriors
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.settings import LOGGER

from ._recursively_test_clade_adherence import recursively_test_clade_adherence
from ._perform_adherence_test_of_child_clades import (
    perform_adherence_test_of_child_clades,
)


def perform_phylogenetic_adherence_test(
    target_sequence: str,
    reference_set: ReferenceSet,
    tree_priors: TreePriors,
) -> Either[c_exc.MappedErrors, bool]:
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

        # TODO: do implement

        # ? --------------------------------------------------------------------
        # ? Collect hierarchical tree
        # ? --------------------------------------------------------------------

        tree_either = reference_set.get_hierarchical_tree()

        if tree_either.is_left:
            return tree_either

        tree: CladeWrapper = tree_either.value

        print(f"Tree: {tree}")

        if tree.is_root() is False:
            return left(
                c_exc.UseCaseError(
                    "Unexpected error. Retrieved tree using "
                    + "`get_hierarchical_tree` is not a rooted tree.",
                    logger=LOGGER,
                )
            )

        # print(tree.get_pretty_tree())

        # ? --------------------------------------------------------------------
        # ? Perform adherence test for outgroup
        # ? --------------------------------------------------------------------

        """ outgroup_adherence_test_either = (
            do_clade_adherence_test_for_single_sequence(
                target_sequence=target_sequence,
                clade_priors=tree_priors.outgroup,
                kmer_indices=reference_set.msa.kmers_indices,
            )
        )

        if outgroup_adherence_test_either.is_left:
            return outgroup_adherence_test_either

        outgroup_adherence_test = outgroup_adherence_test_either.value """

        # ? --------------------------------------------------------------------
        # ? Perform adherence test for the outgroup-sister pairs
        # ? --------------------------------------------------------------------

        first_level_ingroup_clades_either = tree.get_ingroup_clades()

        if first_level_ingroup_clades_either.is_left:
            return first_level_ingroup_clades_either

        first_level_ingroup_clades: list[CladeWrapper] = [
            ingroup
            for ingroup in first_level_ingroup_clades_either.value
            if ingroup.parent == tree_priors.outgroup.parent
        ]

        # ? --------------------------------------------------------------------
        # ? Perform adherence test for the ingroups
        # ? --------------------------------------------------------------------

        for clade in first_level_ingroup_clades:
            print(f"clade: {clade}\n")

            if (children := clade.children) is None:
                continue

            perform_adherence_test_of_child_clades(
                target_sequence=target_sequence,
                clades=[i for i in children if i.is_internal()],
                tree_priors=tree_priors,
                kmer_indices=reference_set.msa.kmers_indices,
                total_length=len(reference_set.labels_map),
            )

            continue

            adherence_response_either = recursively_test_clade_adherence(
                target_sequence=target_sequence,
                clades=[i for i in children if i.is_internal()],
                tree_priors=tree_priors,
                kmer_indices=reference_set.msa.kmers_indices,
                total_length=len(reference_set.labels_map),
            )

            if adherence_response_either.is_left:
                return adherence_response_either

            # print(adherence_response_either.value)

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(True)

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))
