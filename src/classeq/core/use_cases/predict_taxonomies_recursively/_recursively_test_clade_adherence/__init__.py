from typing import Iterator

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


def recursively_test_clade_adherence(
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
                    i for i in tree_priors.ingroups if i.parent == clade.parent
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

            print(clade_priors)

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
                    return c_exc.UseCaseError(
                        "Unexpected error on try to calculate "
                        + f"adherence test for group {name}.",
                        logger=LOGGER,
                    )()

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
            return c_exc.UseCaseError(exc, logger=LOGGER)()

    try:
        LOGGER.debug(f"Received clades: {clades}")

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
        return c_exc.UseCaseError(exc, logger=LOGGER)()
