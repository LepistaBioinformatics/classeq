from copy import deepcopy

import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.domain.dtos.clade import CladeWrapper
from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.priors import PriorGroup, TreePriors
from classeq.settings import LOGGER

from ._do_clade_adherence_test_for_single_sequence import (
    do_clade_adherence_test_for_single_sequence,
)


def perform_adherence_test_of_child_clades(
    target_sequence: str,
    clades: list[CladeWrapper],
    tree_priors: TreePriors,
    kmer_indices: KmersInverseIndices,
    total_length: int,
) -> Either[c_exc.MappedErrors, CladeWrapper]:
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

        contrasting_clades: list[tuple[CladeWrapper, float, float]] = []

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
                contrasting_clades.append((deepcopy(clade), ingroup, sister))

        print(contrasting_clades)

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(True)

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
