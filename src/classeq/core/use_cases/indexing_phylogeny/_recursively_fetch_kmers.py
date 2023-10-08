import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.domain.dtos.clade import ClasseqClade
from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.ordered_tuple import OrderedTuple
from classeq.core.domain.dtos.priors import (
    IngroupPriors,
    IngroupCladePriors,
    SisterCladePriors,
    TreePriors,
)
from classeq.core.use_cases.shared.fetch_kmers_indices import (
    fetch_kmers_indices,
)
from classeq.settings import LOGGER

from ._get_terminal_nodes import get_terminal_nodes


def recursively_fetch_kmers(
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

        ingroups_priors: list[IngroupPriors] = []

        for clade in [i for i in ingroups if i.is_internal()]:
            LOGGER.info(f"Calculating priors for clade: {clade}")

            # ? ----------------------------------------------------------------
            # ? Collect group labels
            # ? ----------------------------------------------------------------

            LOGGER.debug("\tCollecting group labels")

            ingroup: list[ClasseqClade] = get_terminal_nodes(
                target_nodes=[i for i in ingroups if i.parent == clade.id],
                reference_nodes=ingroups,
            )

            if len(ingroup) < min_clade_size:
                LOGGER.debug(
                    "Ingroup clade ineligible for training due do not having "
                    + f"the minimum number of terminals ({min_clade_size}): "
                    + f"{clade.id}"
                )

                continue

            ingroup_labels: list[int] = [i.name for i in ingroup]

            ingroup_ids = [i.id for i in ingroup]

            sister_group: list[ClasseqClade] = [
                terminal
                for terminal in get_terminal_nodes(
                    reference_nodes=ingroups,
                    target_nodes=[
                        i for i in ingroups if i.parent == clade.parent
                    ],
                )
                if terminal.id not in ingroup_ids
            ]

            if len(sister_group) < min_clade_size:
                LOGGER.debug(
                    "Sister group clade ineligible for training due do not "
                    + "reaches the minimum number of terminals "
                    + f"({min_clade_size}): {clade.id}"
                )

                continue

            sister_group_labels: list[int] = [i.name for i in sister_group]

            LOGGER.debug("\tCollection finished")

            # ? ----------------------------------------------------------------
            # ? Estimate specific priors
            # ? ----------------------------------------------------------------

            LOGGER.debug("\tEstimating kmer specific priors")

            LOGGER.debug("\t\tEstimating ingroup")

            if (
                ingroup_priors_either := fetch_kmers_indices(
                    kmer_indices=kmer_indices,
                    sequence_codes=ingroup_labels,
                )
            ).is_left:
                return ingroup_priors_either

            ingroup_specific_priors = IngroupCladePriors(
                labels=OrderedTuple(ingroup_labels),
                kmers=ingroup_priors_either.value,
            )

            LOGGER.debug("\t\tEstimating sister ingroup")

            if (
                sister_group_priors_either := fetch_kmers_indices(
                    kmer_indices=kmer_indices,
                    sequence_codes=sister_group_labels,
                )
            ).is_left:
                return sister_group_priors_either

            sister_group_specific_priors = SisterCladePriors(
                labels=OrderedTuple(sister_group_labels),
                kmers=sister_group_priors_either.value,
            )

            LOGGER.debug("\tEstimation finished")

            ingroups_priors.append(
                IngroupPriors(
                    parent=clade.id,
                    clade_priors=(
                        ingroup_specific_priors,
                        sister_group_specific_priors,
                    ),
                )
            )

        LOGGER.info("Ingroup priors calculated")
        LOGGER.info(f"\t{len(ingroups_priors)} clades eligible for predictions")

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(TreePriors(ingroups=ingroups_priors))

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
