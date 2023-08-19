from uuid import UUID

import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.domain.dtos.clade import ClasseqClade
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
from Bio.Phylo.BaseTree import Tree
from classeq.core.domain.dtos.tree import ClasseqTree
from classeq.settings import LOGGER

from ._estimate_clade_kmer_specific_priors import (
    estimate_clade_kmer_specific_priors,
)
from ._get_terminal_nodes import get_terminal_nodes


def calculate_recursive_priors(
    reference_tree: ClasseqTree,
    root: ClasseqClade,
    outgroups: list[ClasseqClade],
    ingroups: list[ClasseqClade],
    kmer_indices: KmersInverseIndices,
    min_clade_size: int,
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

        if reference_tree.sanitized_tree is None:
            if (
                load_response := reference_tree.load_tree(
                    outgroups=reference_tree.outgroups
                )
            ).is_left:
                return load_response

        sanitized_tree: Tree = reference_tree.sanitized_tree  # type: ignore

        outgroup_clades = sanitized_tree.common_ancestor(
            [i for i in sanitized_tree.get_terminals() if i in outgroups]
        )

        outgroup_path = sanitized_tree.get_path(outgroup_clades)

        if len(outgroup_path) == 0:
            LOGGER.warning("Outgroup is the current tree root")

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

            ingroup: list[ClasseqClade] = get_terminal_nodes(
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

            sister_group: list[ClasseqClade] = [
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
