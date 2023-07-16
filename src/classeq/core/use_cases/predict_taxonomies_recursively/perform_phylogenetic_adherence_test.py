from enum import Enum
from hashlib import md5
from typing import Iterator

import clean_base.exceptions as c_exc
from attrs import define, field
from clean_base.either import Either, left, right

from classeq.core.domain.dtos.clade import CladeWrapper
from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.priors import PriorGroup, TreePriors
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.settings import LOGGER

from .do_clade_adherence_test_for_single_sequence import (
    do_clade_adherence_test_for_single_sequence,
)


class CladeAdherenceResultStatus(Enum):
    """The status of the clade adherence test.

    Description:
        This class defines the status of the clade adherence test. The status
        can be one of the following:
            - `next-iteration`: The test is inconclusive and the test should
                continue to the next iteration.
            - `conclusive-ingroup`: The test is conclusive and the sequence is
                in the ingroup.
            - `conclusive-sister`: The test is conclusive and the sequence is
                in the sister group.
            - `inconclusive`: The test is inconclusive and the sequence is not
                in the ingroup or sister group.
            - `not-possible`: The test is not possible because the clade is a
                leaf node.
            - `unavailable`: The test is not available because the clade is not
                in the tree.

    Attributes:
        NEXT_ITERATION (str): The test is inconclusive and the test should
            continue to the next iteration.
        CONCLUSIVE_INGROUP (str): The test is conclusive and the sequence is
            in the ingroup.
        CONCLUSIVE_SISTER (str): The test is conclusive and the sequence is
            in the sister group.
        INCONCLUSIVE (str): The test is inconclusive and the sequence is not
            in the ingroup or sister group.
        NOT_POSSIBLE (str): The test is not possible because the clade is a
            leaf node.
        UNAVAILABLE (str): The test is not available because the clade is not
            in the tree.

    """

    NEXT_ITERATION = "next-iteration"
    CONCLUSIVE_INGROUP = "conclusive-ingroup"
    CONCLUSIVE_SISTER = "conclusive-sister"
    INCONCLUSIVE = "inconclusive"
    NOT_POSSIBLE = "not-possible"
    UNAVAILABLE = "unavailable"


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

        # print(tree)

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
            print(f"first_level_ingroup_clades clade: {clade}\n")

            if (children := clade.children) is None:
                continue

            adherence_response_either = __recursive_test_clade_adherence(
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


@define
class CladeAdherenceResult:
    clade: CladeWrapper = field()
    parent: set[CladeWrapper] = field()
    status: CladeAdherenceResultStatus = field()
    ingroup_joint_probability: float = field(default=-999)

    # ? ------------------------------------------------------------------------
    # ? Life cycle hook methods
    # ? ------------------------------------------------------------------------

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, CladeAdherenceResult):
            return NotImplemented
        return self.__hash__() == other.__hash__()

    def __ne__(self, other: object) -> bool:
        # Not strictly necessary, but to avoid having both x==y and x!=y True at
        # the same time.
        if not isinstance(other, CladeAdherenceResult):
            return NotImplemented
        return not (self.__hash__() == other.__hash__())

    def __hash__(self) -> int:
        return int(
            md5(
                "-".join(
                    [
                        self.clade.__hash__().__str__(),
                        "".join([i.__str__() for i in self.parent]),
                        self.status.name,
                    ]
                ).encode("utf-8")
            ).hexdigest(),
            16,
        )


def __recursive_test_clade_adherence(
    target_sequence: str,
    clades: list[CladeWrapper],
    tree_priors: TreePriors,
    kmer_indices: KmersInverseIndices,
    total_length: int,
) -> Either[c_exc.MappedErrors, CladeWrapper]:
    """Run clade adherence test for a full tree chain

    Args:
        target_sequence (str): A DNA sequence to be tested.
        clades (list[CladeWrapper]): The clades list to perform tests
        tree_priors (TreePriors): The tree prior object to base tests on.
        kmer_indices (KmersInverseIndices): The kmer indices.

    Returns:
        Either[c_exc.MappedErrors, CladeWrapper]: A single clade if all the
            use-case execution was successfully done.
    """

    def do_single_adherence_test(
        clade: CladeWrapper,
        parent_clades: set[CladeWrapper],
    ) -> Either[c_exc.MappedErrors, tuple[str | None, CladeAdherenceResult]]:
        try:
            try:
                clade_priors = next(
                    i for i in tree_priors.ingroups if i.parent == clade.id
                )
            except StopIteration:
                LOGGER.debug(f"Ignore child: {clade}")
                return right(
                    (
                        "No available priors to test",
                        CladeAdherenceResult(
                            clade=clade,
                            parent=parent_clades,
                            status=CladeAdherenceResultStatus.NOT_POSSIBLE,
                        ),
                    )
                )

            adherence_test_either = do_clade_adherence_test_for_single_sequence(
                target_sequence=target_sequence,
                clade_priors=clade_priors,
                kmer_indices=kmer_indices,
                total_length=total_length,
            )

            if adherence_test_either.is_left:
                return adherence_test_either

            adherence_test = adherence_test_either.value
            ingroup = adherence_test.pop(PriorGroup.INGROUP)
            sister = adherence_test.pop(PriorGroup.SISTER)

            for name, result in [
                ("ingroup", ingroup),
                ("sister", sister),
            ]:
                if result is None:
                    return left(
                        c_exc.UseCaseError(
                            "Unexpected error on try to calculate "
                            + f"adherence test for group {name}.",
                            logger=LOGGER,
                        )
                    )

            # Evaluate if the ingroup probability is greater than others. Case
            # True recursively run adherence test for children clades.
            if ingroup < sister:
                print(clade, ingroup < sister, ingroup, sister)
                print("entrou no ingroup")
                # parent_clades.add(clade)

                if clade.children is None:
                    return right(
                        (
                            "Clade has no children clades to continue",
                            CladeAdherenceResult(
                                clade=clade,
                                parent=parent_clades,
                                status=CladeAdherenceResultStatus.CONCLUSIVE_INGROUP,
                            ),
                        )
                    )

                return right(
                    (
                        None,
                        CladeAdherenceResult(
                            clade=clade,
                            parent=parent_clades,
                            status=CladeAdherenceResultStatus.NEXT_ITERATION,
                            ingroup_joint_probability=ingroup,
                        ),
                    )
                )

            return right(
                (
                    None,
                    CladeAdherenceResult(
                        clade=clade,
                        parent=parent_clades,
                        status=CladeAdherenceResultStatus.CONCLUSIVE_SISTER,
                    ),
                )
            )

        except Exception as exc:
            return left(c_exc.UseCaseError(exc, logger=LOGGER))

    try:
        LOGGER.debug(f"Received clades: {clades}")

        if len(clades) == 0:
            return left(
                c_exc.UseCaseError(
                    "Clades should be an array with at last one element.",
                    logger=LOGGER,
                )
            )

        # ? --------------------------------------------------------------------
        # ? Run adherence test for all clades
        #
        # Adherence test should be executed to all input clades and results
        # stored to further comparisons. Given the possibility to more than one
        # clade would selected all clades should be tested and the output
        # results further evaluated.
        #
        # ? --------------------------------------------------------------------

        def adherence_test_generator(
            clades: list[CladeWrapper],
            parent_clades: set[CladeWrapper],
        ) -> Iterator[tuple[str | None, CladeAdherenceResult]]:
            for clade in clades:
                adherence_response_either = do_single_adherence_test(
                    clade=clade,
                    parent_clades=parent_clades,
                )

                if adherence_response_either.is_left:
                    raise Exception(
                        f"Unexpected error on try to find identity of: {clade}"
                    )

                yield adherence_response_either.value

        """ clades_adherence_responses: list[
            tuple[CladeWrapper, dict[PriorGroup, float]]
        ] = [] """

        status = CladeAdherenceResultStatus.NEXT_ITERATION
        target_clades = clades
        parent_clades: set[CladeWrapper] = set()

        while status == CladeAdherenceResultStatus.NEXT_ITERATION:
            next_clades = {
                result
                for _, result in adherence_test_generator(
                    clades=target_clades,
                    parent_clades=parent_clades,
                )
                if result.status == CladeAdherenceResultStatus.NEXT_ITERATION
            }

            if len(next_clades) == 0:
                print(f"target_clades: {target_clades}")
                print("PASSOU NO BREAK")
                break

            print()
            print("PASSOU DIRETO")
            for next_clade in next_clades:
                print(
                    f"next_clade: ===> {next_clade.ingroup_joint_probability}\t{next_clade}"
                )
            print()

            max_probability_clade = max(
                [clade.ingroup_joint_probability for clade in next_clades]
            )

            unique_target_clades: set[CladeWrapper] = set()
            for group in [
                clade
                for clade in next_clades
                if clade.ingroup_joint_probability == max_probability_clade
            ]:
                for child in [
                    i for i in group.clade.children if i.is_internal()
                ]:
                    unique_target_clades.add(child)

            print()
            for c in unique_target_clades:
                print(f"unique_target_clade: {c}")
            print()

            target_clades = list(set(unique_target_clades))

            if len(target_clades) == 0:
                continue

            print()
            print("PASSOU AQUI")

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(True)

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))
