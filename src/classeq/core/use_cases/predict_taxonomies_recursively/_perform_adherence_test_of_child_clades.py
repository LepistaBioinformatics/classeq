from copy import deepcopy
from enum import Enum

from hashlib import md5
from attrs import define, field
import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.domain.dtos.clade import CladeWrapper
from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.priors import PriorGroup, TreePriors
from classeq.settings import LOGGER

from ._do_clade_adherence_test_for_single_sequence import (
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


@define
class CladeAdherenceResult:
    clade: CladeWrapper = field()
    ingroup_joint_probability: float = field(default=-999)
    sister_joint_probability: float = field(default=-999)

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
                        str(self.ingroup_joint_probability),
                        str(self.sister_joint_probability),
                    ]
                ).encode("utf-8")
            ).hexdigest(),
            16,
        )


def perform_adherence_test_of_child_clades(
    target_sequence: str,
    clades: list[CladeWrapper],
    tree_priors: TreePriors,
    kmer_indices: KmersInverseIndices,
    total_length: int,
) -> Either[
    c_exc.MappedErrors,
    tuple[CladeAdherenceResult | None, CladeAdherenceResultStatus],
]:
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

        # ! fazer o teste de aderência para todas as clades filhos e retornar o com maior pontuação

        contrasting_clades: set[CladeAdherenceResult] = set()

        for clade in clades:
            try:
                clade_priors = next(
                    i for i in tree_priors.ingroups if i.parent == clade.id
                )
            except StopIteration:
                LOGGER.debug(f"Ignore child: {clade}")
                raise NotImplementedError()

            adherence_test_either = do_clade_adherence_test_for_single_sequence(
                target_sequence=target_sequence,
                clade_priors=clade_priors,
                kmer_indices=kmer_indices,
                total_length=total_length,
            )

            if adherence_test_either.is_left:
                return adherence_test_either

            adherence_test = adherence_test_either.value

            if (ingroup := adherence_test.pop(PriorGroup.INGROUP)) is None:
                return c_exc.UseCaseError(
                    "Unexpected error on try to calculate "
                    + "adherence test for ingroup.",
                    logger=LOGGER,
                )()

            if (sister := adherence_test.pop(PriorGroup.SISTER)) is None:
                return c_exc.UseCaseError(
                    "Unexpected error on try to calculate "
                    + "adherence test for sister.",
                    logger=LOGGER,
                )()

            if ingroup < sister:
                contrasting_clades.add(
                    CladeAdherenceResult(
                        clade=deepcopy(clade),
                        ingroup_joint_probability=ingroup,
                        sister_joint_probability=sister,
                    )
                )

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        if len(contrasting_clades) == 0:
            return (
                None,
                CladeAdherenceResultStatus.NOT_POSSIBLE,
            )

        if len(contrasting_clades) == 1:
            clade_binding: CladeAdherenceResult

            if (clade_binding := contrasting_clades.pop()) is None:
                return c_exc.UseCaseError(
                    "Unexpected error on try to get clade from set.",
                    logger=LOGGER,
                )()

            if clade_binding.clade.children is not None:
                return (
                    clade_binding,
                    CladeAdherenceResultStatus.NEXT_ITERATION,
                )

            return (
                clade_binding,
                CladeAdherenceResultStatus.CONCLUSIVE_INGROUP
                if clade_binding.ingroup_joint_probability
                < clade_binding.sister_joint_probability
                else CladeAdherenceResultStatus.CONCLUSIVE_SISTER,
            )

        if len(contrasting_clades) > 1:
            return (
                None,
                CladeAdherenceResultStatus.INCONCLUSIVE,
            )

        return right(True)

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
