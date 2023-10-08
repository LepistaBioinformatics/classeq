import clean_base.exceptions as c_exc
from clean_base.either import Either

from classeq.core.domain.dtos.priors import TreePriors
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.settings import LOGGER

from ._recursively_fetch_kmers import recursively_fetch_kmers


def fetch_clade_specific_kmers(
    references: ReferenceSet,
    min_clade_size: int,
) -> Either[c_exc.MappedErrors, TreePriors]:
    """Estimate clade specific priors.

    Description:
        This function estimates the clade specific priors for each clade in the
        phylogenetic tree. The function is recursive and it calculates the
        probabilities for each clade in the tree.

    Args:
        references (ReferenceSet): Reference set with the phylogenetic tree and
            the multiple sequence alignment.

    Returns:
        Either[c_exc.MappedErrors, TreePriors]: Either with the clade specific
            priors or an error.

    Raises:
        c_exc.UseCaseError: If the argument `references` is not a
            `ReferenceSet` instance.
        c_exc.UseCaseError: If the linear tree could not be generated or if the
            root node is not present in the tree.

    """

    try:
        # ? --------------------------------------------------------------------
        # ? Validate entries
        # ? --------------------------------------------------------------------

        if not isinstance(references, ReferenceSet):
            return c_exc.UseCaseError(
                f"Argument `references` should be a `{ReferenceSet}` instance.",
                exp=True,
                logger=LOGGER,
            )()

        # ? --------------------------------------------------------------------
        # ? Get linear tree
        #
        # Different from default `Bio.Phylo.BaseTree.Tree` from 'BioPython's'
        # library, classeq linear trees contain indication for the parent node
        # of each node. This link allows to perform recursive filtration along
        # the tree. This is usual over likelihood calculations.
        #
        # ? --------------------------------------------------------------------

        if references.linear_tree is None:
            linear_tree_either = references.build_linear_tree()

            if linear_tree_either.is_left:
                return linear_tree_either

        if references.linear_tree is None:
            return c_exc.UseCaseError(
                "Could not generate linear tree.",
                logger=LOGGER,
            )()

        # ? --------------------------------------------------------------------
        # ? Recursive calculate probabilities
        # ? --------------------------------------------------------------------

        return recursively_fetch_kmers(
            ingroups=references.linear_tree,
            kmer_indices=references.msa.kmers_indices,
            min_clade_size=min_clade_size,
        )

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
