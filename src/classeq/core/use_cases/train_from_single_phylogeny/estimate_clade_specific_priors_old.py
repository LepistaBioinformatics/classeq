from collections import defaultdict
from functools import reduce
from typing import DefaultDict, Iterator
from classeq.core.domain.dtos.clade import CladeWrapper
from classeq.core.domain.dtos.kmer_inverse_index import (
    KmerIndex,
    KmersInverseIndices,
)
import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import DEFAULT_KMER_SIZE, LOGGER


def estimate_clade_specific_priors(
    references: ReferenceSet,
) -> Either[c_exc.MappedErrors, bool]:
    """# ! DEPRECATED"""

    try:
        # ? --------------------------------------------------------------------
        # ? Validate entries
        # ? --------------------------------------------------------------------

        # TODO

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

        linear_tree: set[CladeWrapper] = linear_tree_either.value

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

        __calculate_recursive_priors(
            root=root,
            outgroups=outgroup_nodes,
            other_nodes=[
                i
                for i in linear_tree
                if i.id not in [root.id, *[o.id for o in outgroup_nodes]]
            ],
            label_map=references.labels_map,
            kmer_indices=references.msa.kmers_indices,
        )

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(True)

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))


def __calculate_recursive_priors(
    root: CladeWrapper,
    outgroups: list[CladeWrapper],
    other_nodes: list[CladeWrapper],
    label_map: DefaultDict[str, int],
    kmer_indices: KmersInverseIndices,
) -> Either[c_exc.MappedErrors, float]:
    def __calculate_priors_for_clade(
        sequence_codes: list[int],
    ) -> Either[c_exc.MappedErrors, float]:
        # ? --------------------------------------------------------------------
        # ? Estimate outgroup/sister pair kmer priors
        # ? --------------------------------------------------------------------

        kmer_priors_either = __estimate_clade_kmer_specific_priors(
            kmer_indices=kmer_indices,
            sequence_codes=sequence_codes,
            corpus_size=(1 + len(outgroups) + len(other_nodes)),
        )

        if kmer_priors_either.is_left:
            return kmer_priors_either

        kmer_priors: DefaultDict[str, float] = kmer_priors_either.value

        # ? --------------------------------------------------------------------
        # ? Calculate joint probability units
        # ? --------------------------------------------------------------------

        joint_probability_units: list[float] = []

        for kmer in target_sequence_kmers:
            if (prior := kmer_priors.get(kmer)) is None:
                continue

            kmer_index = kmer_indices.index_of(kmer)

            if kmer_index is None:
                return left(
                    c_exc.UseCaseError(
                        "Unexpected error on calculate outgroup joint "
                        + f"probability. `{kmer}` is not indexed.",
                        logger=LOGGER,
                    )
                )

            current_index = [
                i
                for i in sequence_codes
                if kmer_indices.indices[kmer_index].contains(i)
            ]

            joint_probability_units.append(
                __calculate_probability_of_group_contains_kmer(
                    prior=prior,
                    sequences_with_kmer=len(current_index),
                    total_sequences=len(sequence_codes),
                )
            )

        # ? --------------------------------------------------------------------
        # ? Built the joint probability product
        # ? --------------------------------------------------------------------

        print(f"joint_probability_units: {joint_probability_units}")

        return right(reduce(lambda i, j: i * j, joint_probability_units))

    try:
        # ? --------------------------------------------------------------------
        # ? Target sequence
        # ? --------------------------------------------------------------------

        target_sequence_kmers = [
            kmer
            for kmer in KmersInverseIndices.generate_kmers(
                # Col_cuscutae_CSL_473
                # dna_sequence="CCTTCATTGAGACCAAGTACGCTGTGAGTATCACCCCCACTTTACCCCTCCATAATGATATCACGTCTGCTACAATAACACCAGCTTCATCGGTAACCACGGGAAAAGAGTCAGAGCTAGTACTCTCGACTCTTTGGACCCAAGGTTTCGATTGGGCTCGTTGTTGTAATGATACGACGTGACACAATCATGCAGAAACAGCCCAAACAAAATTTGCTGACAGACAATCATCACAGGCCTACATGCTCAAGTAC",
                # Col_orchidophilum_IMI_309357
                # dna_sequence="CCTTCATTGAGACCAAGTACGCTGTGAGTATCACCCCACTTTACCCCTCCATGATGATATCACATCTGTCACGACAATACCAGCCTCATCGGCCACTGGGAAAGAAATGAGCTAGCACTCTCGATCCTGTGACCCAGGATACTGAAGCGGCTCGTCCCAATGGCATGATGTGACTAGGTCACGAAGAAATAGTTGGGACAACATTTGCTGACAGACCACTACCACAGGCCTACATGCTCAAGTAC",
                # Col_acutatum_CBS_110735
                dna_sequence="CCTTCATTGAGACCAAGTACGCTGTGAGTATCACCCCACTTTACCCCTCCATCATGATATCACGTCTGCCACGATAACACCAGCTTCGTCGGTACCCACGGGAAAAGAGTCAGAGCTAGCGCTCTCGACTCTTTTGCCCCGAGGTTTCGATTGGGCTCGTTGTAATGATGCGACGTGATACAACCATGCAGAAACAGCCGAGACAAAATTTGCTGACAGACAATCATCACAGGCCTACATGCTCAAGTAC",
                k_size=DEFAULT_KMER_SIZE,
            )
        ]

        # ? --------------------------------------------------------------------
        # ? Calculate outgroup joint probabilities
        # ? --------------------------------------------------------------------

        outgroup_labels: list[int] = []

        for node in outgroups:
            if (label := label_map.get(node.name)) is None:
                return left(
                    c_exc.UseCaseError(
                        f"Invalid labels map. `{node}` is not included.",
                        logger=LOGGER,
                    )
                )

            outgroup_labels.append(label)

        outgroup_joint_probability_either = __calculate_priors_for_clade(
            sequence_codes=outgroup_labels
        )

        if outgroup_joint_probability_either.is_left:
            return outgroup_joint_probability_either

        outgroup_joint_probability: float = (
            outgroup_joint_probability_either.value
        )

        # ? --------------------------------------------------------------------
        # ? Calculate outgroup-sister nodes likelihood
        # ? --------------------------------------------------------------------

        outgroup_sisters_labels: list[int] = []

        for i in __get_terminal_nodes(
            target_nodes=[i for i in other_nodes if i.parent == root.id],
            reference_nodes=other_nodes,
        ):
            if (label := label_map.get(i.name)) is None:
                return left(
                    c_exc.UseCaseError(
                        f"Invalid labels map. `{node}` is not included.",
                        logger=LOGGER,
                    )
                )

            outgroup_sisters_labels.append(label)

        sister_joint_probability_either = __calculate_priors_for_clade(
            sequence_codes=outgroup_sisters_labels
        )

        if sister_joint_probability_either.is_left:
            return sister_joint_probability_either

        sister_joint_probability: float = sister_joint_probability_either.value

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        outgroup_joint_probability = -outgroup_joint_probability
        sister_joint_probability = -sister_joint_probability

        print(f"outgroup_joint_probability: {outgroup_joint_probability}")
        print(f"sister_joint_probability: {sister_joint_probability}")
        print(outgroup_joint_probability > sister_joint_probability)
        print(sister_joint_probability > outgroup_joint_probability)

        return right(True)

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))


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


def __estimate_clade_kmer_specific_priors(
    kmer_indices: KmersInverseIndices,
    sequence_codes: list[int],
    corpus_size: int,
) -> Either[c_exc.MappedErrors, DefaultDict[str, float]]:
    try:
        kmers_priors_for_clade: DefaultDict[str, float] = defaultdict()
        target_indices: set[KmerIndex] = set()

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
            kmers_priors_for_clade[kmer_index.kmer] = prior

        return right(kmers_priors_for_clade)

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))


def __get_terminal_nodes(
    target_nodes: list[CladeWrapper],
    reference_nodes: list[CladeWrapper],
) -> list[CladeWrapper]:
    def __get_children(
        target_nodes: list[CladeWrapper],
    ) -> Iterator[CladeWrapper]:
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
