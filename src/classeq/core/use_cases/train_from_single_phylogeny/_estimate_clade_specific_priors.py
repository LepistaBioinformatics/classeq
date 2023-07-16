from uuid import UUID

import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.domain.dtos.clade import CladeWrapper
from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.priors import (
    IngroupCladePriors,
    IngroupLabeledPriors,
    LabeledPriors,
    OutgroupCladePriors,
    OutgroupLabeledPriors,
    PriorGroup,
    SisterGroupLabeledPriors,
    TreePriors,
)
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.settings import LOGGER

from ._estimate_clade_kmer_specific_priors import (
    estimate_clade_kmer_specific_priors,
)
from ._get_terminal_nodes import get_terminal_nodes


def estimate_clade_specific_priors(
    references: ReferenceSet,
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
        c_exc.InvalidArgumentError: If the argument `references` is not a
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
        # ? Extract root node
        #
        # Recursive analysis starts from the root node. Them the root should be
        # identified before the analysis initialization.
        #
        # ? --------------------------------------------------------------------

        try:
            root = next(i for i in iter(references.linear_tree) if i.is_root())
        except StopIteration:
            return c_exc.UseCaseError(
                "Root node not present in linear tree.",
                logger=LOGGER,
            )()

        # ? --------------------------------------------------------------------
        # ? Extract outgroup nodes
        #
        # Outgroup nodes should be included over each training step as a random
        # noise.
        #
        # ? --------------------------------------------------------------------

        outgroup_nodes = [i for i in references.linear_tree if i.is_outgroup()]
        expected_outgroups = [
            name
            for i in references.tree.outgroups
            if (name := references.labels_map.get(i)) is not None
        ]

        if not all([i.name in expected_outgroups for i in outgroup_nodes]):
            return c_exc.UseCaseError(
                "Not all outgroups are present at the phylogenetic tree. "
                + f"Expected nodes: {', '.join(references.tree.outgroups)}",
                logger=LOGGER,
            )()

        # ? --------------------------------------------------------------------
        # ? Recursive calculate probabilities
        # ? --------------------------------------------------------------------

        return __calculate_recursive_priors(
            root=root,
            outgroups=outgroup_nodes,
            ingroups=[
                i
                for i in references.linear_tree
                if i.id not in [root.id, *[o.id for o in outgroup_nodes]]
            ],
            kmer_indices=references.msa.kmers_indices,
        )

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()


def __calculate_recursive_priors(
    root: CladeWrapper,
    outgroups: list[CladeWrapper],
    ingroups: list[CladeWrapper],
    kmer_indices: KmersInverseIndices,
    min_clade_size: int = 1,
) -> Either[c_exc.MappedErrors, TreePriors]:
    """Recursive calculate probabilities for each clade.

    Description:
        This function is the core of the recursive calculation of the
        probabilities for each clade. The function is recursive and it
        calculates the probabilities for each clade in the tree.

    Args:
        root (CladeWrapper): Root node of the tree.
        outgroups (list[CladeWrapper]): Outgroup nodes of the tree.
        ingroups (list[CladeWrapper]): Ingroup nodes of the tree.
        kmer_indices (KmersInverseIndices): Kmer indices of the MSA.
        min_clade_size (int, optional): Minimum size of the clade to be
            considered. Defaults to 1.

    Returns:
        Either[c_exc.MappedErrors, TreePriors]: Either a `TreePriors` instance
            or a `classeq.core.domain.utils.exceptions.MappedErrors` instance.

    Raises:
        c_exc.UseCaseError: If the outgroups are not valid.

    """

    try:
        # ? --------------------------------------------------------------------
        # ? Calculate outgroup clade priors
        #
        # Outgroup specific priors should be calculated separately to test
        # specific hypothesis during the prediction steps.
        #
        # ? --------------------------------------------------------------------

        outgroup_labels: list[int] = []
        outgroups_parents = {o.parent for o in outgroups}

        if len(outgroups_parents) != 1:
            return c_exc.UseCaseError(
                "Invalid outgroups. Outgroups has more than one parent node"
                + f": {outgroups_parents}",
                logger=LOGGER,
            )()

        outgroup_parent: UUID = outgroups_parents.pop()

        if root.id != outgroup_parent:
            return c_exc.UseCaseError(
                "Invalid outgroups. Outgroups should share the the root as "
                + f"parent node. Expected {root.id}, found {outgroup_parent}",
                logger=LOGGER,
            )()

        outgroup_priors_either = estimate_clade_kmer_specific_priors(
            kmer_indices=kmer_indices,
            sequence_codes=[i.name for i in outgroups],
            corpus_size=(1 + len(outgroups) + len(ingroups)),
        )

        if outgroup_priors_either.is_left:
            return outgroup_priors_either

        outgroup_priors = OutgroupCladePriors(
            parent=outgroup_parent,
            priors=OutgroupLabeledPriors(
                labels=tuple(sorted(outgroup_labels)),
                priors=outgroup_priors_either.value,
            ),
        )

        # ? --------------------------------------------------------------------
        # ? Calculate ingroups clades priors
        #
        # Ingroup specific priors including a concatenation of three subparts:
        #
        # Ingroup: including only kmer priors of the target clade.
        #
        # Sister group: including kmer priors of the all sister clades of the
        # ingroup clade.
        #
        # Random group: including a sample of terminals not present on ingroup
        # or sister group. Random elements should be used to estimate the
        # probability of the prediction sequence to does not belongs to the
        # current super-clade.
        #
        # ? --------------------------------------------------------------------

        ingroups_priors: list[IngroupCladePriors] = []

        for clade in [i for i in ingroups if i.is_internal()]:
            # ? ----------------------------------------------------------------
            # ? Extract ingroup terminals
            # ? ----------------------------------------------------------------

            ingroup: list[CladeWrapper] = get_terminal_nodes(
                target_nodes=[i for i in ingroups if i.parent == clade.id],
                reference_nodes=ingroups,
            )

            if len(ingroup) < min_clade_size:
                LOGGER.warning(
                    "Ingroup clade ineligible for training due to not having "
                    + f"the minimum number of terminals ({min_clade_size}): "
                    + f"{clade.id}"
                )

                continue

            # ? ----------------------------------------------------------------
            # ? Extract sister group terminals
            # ? ----------------------------------------------------------------

            sister_group: list[CladeWrapper] = [
                terminal
                for terminal in get_terminal_nodes(
                    target_nodes=[
                        i for i in ingroups if i.parent == clade.parent
                    ],
                    reference_nodes=ingroups,
                )
                if terminal.id not in [i.id for i in ingroup]
            ]

            if len(sister_group) < min_clade_size:
                LOGGER.warning(
                    "Sister group clade ineligible for training due to not "
                    + "reaches the minimum number of terminals "
                    + f"({min_clade_size}): {clade.id}"
                )

                continue

            # ? ----------------------------------------------------------------
            # ? Estimate group specific priors
            # ? ----------------------------------------------------------------

            priors_values: list[LabeledPriors] = []

            for receiver, group in [
                (IngroupLabeledPriors, ingroup),
                (SisterGroupLabeledPriors, sister_group),
            ]:
                group_labels: list[int] = [i.name for i in group]

                group_priors_either = estimate_clade_kmer_specific_priors(
                    kmer_indices=kmer_indices,
                    sequence_codes=group_labels,
                    # corpus_size=(1 + len(outgroups) + len(ingroups)),
                    corpus_size=(1 + len(ingroup) + len(sister_group)),
                )

                if group_priors_either.is_left:
                    return group_priors_either

                priors_values.append(
                    receiver(
                        labels=tuple(sorted(group_labels)),
                        priors=group_priors_either.value,
                    )
                )

            ingroup_priors = next(
                i for i in priors_values if i.group == PriorGroup.INGROUP
            )

            sister_priors = next(
                i for i in priors_values if i.group == PriorGroup.SISTER
            )

            ingroups_priors.append(
                IngroupCladePriors(
                    parent=clade.id,
                    priors=(
                        ingroup_priors,
                        sister_priors,
                    ),
                )
            )

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(
            TreePriors(
                outgroup=outgroup_priors,
                ingroups=ingroups_priors,
            )
        )

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
