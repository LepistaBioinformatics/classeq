from functools import reduce
from typing import List, Tuple

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.priors import (
    IngroupCladePriors,
    LabeledPriors,
    PriorGroup,
)
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import DEFAULT_KMER_SIZE, LOGGER


def do_clade_adherence_test_for_single_sequence(
    target_sequence: str,
    clade_priors: IngroupCladePriors,
    kmer_indices: KmersInverseIndices,
) -> Either[bool, c_exc.MappedErrors]:
    try:
        # ? --------------------------------------------------------------------
        # ? Validate entries
        # ? --------------------------------------------------------------------

        # TODO

        # ? --------------------------------------------------------------------
        # ? Get sequence kmers
        # ? --------------------------------------------------------------------

        target_sequence_kmers = [
            kmer
            for kmer in KmersInverseIndices.generate_kmers(
                dna_sequence=target_sequence,
                k_size=DEFAULT_KMER_SIZE,
            )
        ]

        # ? --------------------------------------------------------------------
        # ? Calculate joint probability units
        # ? --------------------------------------------------------------------

        clade_adherence_stats = List[Tuple[PriorGroup, float]]

        for group in clade_priors.priors:
            group_adherence_either = __calculate_clade_adherence(
                priors=group,
                target_kmers=target_sequence_kmers,
                kmer_indices=kmer_indices,
            )

            if group_adherence_either.is_left:
                return group_adherence_either

            clade_adherence_stats.append(
                group.group,
                group_adherence_either.value,
            )

        print(f"clade_adherence_stats: {clade_adherence_stats}")

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(True)

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))


def __calculate_clade_adherence(
    priors: LabeledPriors,
    target_kmers: List[str],
    kmer_indices: KmersInverseIndices,
) -> Either[bool, c_exc.MappedErrors]:
    try:
        joint_probability_units: List[float] = []

        for kmer in target_kmers:
            if (prior := priors.get(kmer)) is None:
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
                for i in priors.labels
                if kmer_indices.indices[kmer_index].contains(i)
            ]

            joint_probability_units.append(
                __calculate_probability_of_group_contains_kmer(
                    prior=prior,
                    sequences_with_kmer=len(current_index),
                    total_sequences=len(priors.labels),
                )
            )

        print(f"joint_probability_units: {joint_probability_units}")

        return right(reduce(lambda i, j: i * j, joint_probability_units))

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
