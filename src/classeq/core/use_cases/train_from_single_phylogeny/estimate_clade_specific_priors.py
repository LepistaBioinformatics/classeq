from collections import defaultdict
from random import sample
from typing import DefaultDict, Dict, Iterator, List, Set
from uuid import UUID

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.clade import CladeWrapper
from classeq.core.domain.dtos.kmer_inverse_index import (
    KmerIndex,
    KmersInverseIndices,
)
from classeq.core.domain.dtos.priors import (
    CladePriors,
    IngroupCladePriors,
    TreePriors,
)
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER


def estimate_clade_specific_priors(
    references: ReferenceSet,
) -> Either[TreePriors, c_exc.MappedErrors]:
    try:
        # ? --------------------------------------------------------------------
        # ? Validate entries
        # ? --------------------------------------------------------------------

        if not isinstance(references, ReferenceSet):
            return left(
                c_exc.InvalidArgumentError(
                    f"Argument `references` should be a `{ReferenceSet}` instance.",
                    exp=True,
                    logger=LOGGER,
                )
            )

        # ? --------------------------------------------------------------------
        # ? Get linear tree
        #
        # Different from default `Bio.Phylo.BaseTree.Tree` from 'BioPython's'
        # library, classeq linear trees contain indication for the parent node
        # of each node. This link allows to perform recursive filtration along
        # the tree. This is usual over likelihood calculations.
        #
        # ? --------------------------------------------------------------------

        linear_tree_either = references.get_linear_tree()

        if linear_tree_either.is_left:
            return linear_tree_either

        linear_tree: Set[CladeWrapper] = linear_tree_either.value

        # ? --------------------------------------------------------------------
        # ? Extract root node
        #
        # Recursive analysis starts from the root node. Them the root should be
        # identified before the analysis initialization.
        #
        # ? --------------------------------------------------------------------

        try:
            root = next(i for i in iter(linear_tree) if i.is_root())
        except StopIteration:
            return left(
                c_exc.UseCaseError(
                    "Root node not present in linear tree.",
                    logger=LOGGER,
                )
            )

        # ? --------------------------------------------------------------------
        # ? Extract outgroup nodes
        #
        # Outgroup nodes should be included over each training step as a random
        # noise.
        #
        # ? --------------------------------------------------------------------

        outgroup_nodes = [i for i in linear_tree if i.is_outgroup()]

        if not all(
            [i.name in references.tree.outgroups for i in outgroup_nodes]
        ):
            return left(
                c_exc.UseCaseError(
                    "Not all outgroups are present at the phylogenetic tree. "
                    + f"Expected nodes: {', '.join(references.tree.outgroups)}",
                    logger=LOGGER,
                )
            )

        # ? --------------------------------------------------------------------
        # ? Recursive calculate probabilities
        # ? --------------------------------------------------------------------

        return __calculate_recursive_priors(
            root=root,
            outgroups=outgroup_nodes,
            ingroups=[
                i
                for i in linear_tree
                if i.id not in [root.id, *[o.id for o in outgroup_nodes]]
            ],
            label_map=references.labels_map,
            kmer_indices=references.msa.kmers_indices,
        )

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))


def __calculate_recursive_priors(
    root: CladeWrapper,
    outgroups: List[CladeWrapper],
    ingroups: List[CladeWrapper],
    label_map: DefaultDict[str, int],
    kmer_indices: KmersInverseIndices,
    min_clade_size: int = 5,
) -> Either[TreePriors, c_exc.MappedErrors]:
    try:
        # ? --------------------------------------------------------------------
        # ? Calculate outgroup clade priors
        #
        # Outgroup specific priors should be calculated separately to test
        # specific hypothesis during the prediction steps.
        #
        # ? --------------------------------------------------------------------

        outgroup_labels: List[int] = []
        outgroups_parents = {o.parent for o in outgroups}

        if len(outgroups_parents) != 1:
            return left(
                c_exc.UseCaseError(
                    "Invalid outgroups. Outgroups has more than one parent node"
                    + f": {outgroups_parents}",
                    logger=LOGGER,
                )
            )

        outgroup_parent: UUID = outgroups_parents.pop()

        if root.id != outgroup_parent:
            return left(
                c_exc.UseCaseError(
                    "Invalid outgroups. Outgroups should share the the root as "
                    + f"parent node. Expected {root.id}, found {outgroup_parent}",
                    logger=LOGGER,
                )
            )

        for node in outgroups:
            if (label := label_map.get(node.name)) is None:
                return left(
                    c_exc.UseCaseError(
                        f"Invalid labels map. `{node}` is not included.",
                        logger=LOGGER,
                    )
                )

            outgroup_labels.append(label)

        outgroup_priors_either = __estimate_clade_kmer_specific_priors(
            kmer_indices=kmer_indices,
            sequence_codes=outgroup_labels,
            corpus_size=(1 + len(outgroups) + len(ingroups)),
        )

        if outgroup_priors_either.is_left:
            return outgroup_priors_either

        outgroup_priors = CladePriors(
            parent=outgroup_parent,
            priors=outgroup_priors_either.value,
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

        ingroups_priors: List[IngroupCladePriors] = []

        for clade in [i for i in ingroups if i.is_internal()]:
            # ? ----------------------------------------------------------------
            # ? Extract ingroup terminals
            # ? ----------------------------------------------------------------

            ingroup: List[CladeWrapper] = __get_terminal_nodes(
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

            sister_group: List[CladeWrapper] = [
                terminal
                for terminal in __get_terminal_nodes(
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
                    + "having the minimum number of terminals "
                    + f"({min_clade_size}): {clade.id}"
                )

                continue

            # ? ----------------------------------------------------------------
            # ? Extract noise group terminals
            # ? ----------------------------------------------------------------

            remaining = [*ingroup, *sister_group]
            noise_group: List[CladeWrapper] = [
                terminal
                for terminal in __get_terminal_nodes(
                    target_nodes=[root],
                    reference_nodes=[*outgroups, *ingroups],
                )
                if terminal.id not in [i.id for i in remaining]
            ]

            noise_group = sample(
                noise_group,
                len(remaining)
                if len(noise_group) > len(remaining)
                else len(noise_group),
            )

            if len(noise_group) < min_clade_size:
                LOGGER.warning(
                    "Noise group clade ineligible for training due to not "
                    + "having the minimum number of terminals "
                    + f"({min_clade_size}): {clade.id}"
                )

                continue

            # ? ----------------------------------------------------------------
            # ? Estimate group specific priors
            # ? ----------------------------------------------------------------

            priors_values: Dict[str, float | None] = {
                "ingroup_priors": None,
                "sister_group_priors": None,
                "noise_group_priors": None,
            }

            for receiver, group in [
                ("ingroup_priors", ingroup),
                ("sister_group_priors", sister_group),
                ("noise_group_priors", noise_group),
            ]:
                group_labels: List[int] = []

                for item in group:
                    if (label := label_map.get(item.name)) is None:
                        return left(
                            c_exc.UseCaseError(
                                f"Invalid labels map. `{item.name}` is not "
                                + "included.",
                                logger=LOGGER,
                            )
                        )

                    group_labels.append(label)

                group_priors_either = __estimate_clade_kmer_specific_priors(
                    kmer_indices=kmer_indices,
                    sequence_codes=group_labels,
                    corpus_size=(1 + len(outgroups) + len(ingroups)),
                )

                if group_priors_either.is_left:
                    return group_priors_either

                priors_values[receiver] = group_priors_either.value

            for key, value in priors_values.items():
                if value is None:
                    return left(
                        c_exc.UseCaseError(
                            f"Unexpected error on calculate priors of `{key}` "
                            + f"for clade {clade}.",
                            logger=LOGGER,
                        )
                    )

            ingroups_priors.append(
                IngroupCladePriors(
                    parent=clade.parent,
                    ingroup_priors=priors_values.get("ingroup_priors"),
                    sister_group_priors=priors_values.get(
                        "sister_group_priors"
                    ),
                    noise_group_priors=priors_values.get("noise_group_priors"),
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
        return left(c_exc.UseCaseError(exc, logger=LOGGER))


def __estimate_clade_kmer_specific_priors(
    kmer_indices: KmersInverseIndices,
    sequence_codes: List[int],
    corpus_size: int,
) -> Either[DefaultDict[str, float], c_exc.MappedErrors]:
    try:
        kmers_priors_for_clade: DefaultDict[str, float] = defaultdict()
        target_indices: Set[KmerIndex] = set()

        for index in kmer_indices.indices:
            for code in sequence_codes:
                if index.contains(code):
                    target_indices.add(index)
                    continue

        for kmer_index in iter(target_indices):
            n = [i for i in kmer_index.records if i in sequence_codes]

            if len(n) == 0:
                raise

            prior = (len(n) + 0.5) / (corpus_size + 1)
            kmers_priors_for_clade[kmer_index.kmer] = round(prior, 8)

        return right(kmers_priors_for_clade)

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))


def __get_terminal_nodes(
    target_nodes: List[CladeWrapper],
    reference_nodes: List[CladeWrapper],
) -> List[CladeWrapper]:
    """Recursively get nodes of type NodeType.TERMINAL.

    Args:
        target_nodes (List[CladeWrapper]): A list of clades to find terminals
            from.
        reference_nodes (List[CladeWrapper]): A list of the remaining nodes to
            find complementary terminal nodes.

    Returns:
        List[CladeWrapper]: A list containing only NodeType.TERMINAL.
    """

    def __get_children(
        target_nodes: List[CladeWrapper],
    ) -> Iterator[CladeWrapper]:
        """Performs a recursive search for terminal nodes.

        Args:
            target_nodes (List[CladeWrapper]): A list of the clades to search
                terminals from.

        Yields:
            Iterator[CladeWrapper]: A single clade of the type
                NodeType.TERMINAL.
        """

        node: CladeWrapper
        for node in target_nodes:
            if node.is_terminal():
                yield node

            yield from __get_terminal_nodes(
                target_nodes=[
                    i for i in reference_nodes if i.parent == node.id
                ],
                reference_nodes=[i for i in reference_nodes if i.id != node.id],
            )

    return [i for i in __get_children(target_nodes)]


def __calculate_probability_of_group_contains_kmer(
    prior: float,
    sequences_with_kmer: int,
    total_sequences: int,
) -> float:
    if sequences_with_kmer > total_sequences:
        raise Exception(
            "`total_sequences` could not be greater than `sequences_with_kmer`"
        )

    # ? ------------------------------------------------------------------------
    # ? Calculate
    #
    # The original formulation for the probability that a group contains a kmer
    # is:
    #
    # P(wi|G) = [ m(wi) + Pi ] / ( M + 1 )
    #
    # m(wi) = be the number of the group sequences containing word wi (argument
    # `sequences_with_kmer` of these function).
    #
    # Pi = The prior probability for the kmer in the overall training dataset
    # (argument `prior` of these function).
    #
    # M = The dataset size for the target training group (argument
    # `total_sequences` of these function).
    #
    # ? ------------------------------------------------------------------------

    return (sequences_with_kmer + prior) / (total_sequences + 1)
