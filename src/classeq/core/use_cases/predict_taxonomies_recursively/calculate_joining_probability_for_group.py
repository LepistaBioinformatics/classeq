from typing import List, Set

from attr import define, field

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER


@define(kw_only=True)
class __KmerStats:
    group_length: int = field()
    records_on_group: int = field()

    def get_percentage(self) -> float:
        return (self.records_on_group * 100) / self.group_length


@define(kw_only=True)
class KmerBayesStats:
    kmer: str
    true_positives_post_prob: float


def calculate_joining_probability_for_group(
    prediction_target: List[str],
    positive_reference: List[List[str]],
    negative_reference: List[List[str]],
) -> Either[List[KmerBayesStats], c_exc.MappedErrors]:
    try:
        # ? --------------------------------------------------------------------
        # ? Validate args
        # ? --------------------------------------------------------------------

        if not isinstance(prediction_target, list):
            return left(
                c_exc.InvalidArgumentError(
                    f"Prediction target is not a list: {prediction_target}"
                )
            )

        for name, group in [
            ("positive_reference", positive_reference),
            ("negative_reference", negative_reference),
        ]:
            if not isinstance(group, list):
                return left(
                    c_exc.InvalidArgumentError(f"{name} is not a list: {group}")
                )

            for sub_item in group:
                if not isinstance(sub_item, list):
                    return left(
                        c_exc.InvalidArgumentError(
                            f"An element of {name} is not a list: {sub_item}"
                        )
                    )

        # ? --------------------------------------------------------------------
        # ? Calculate joining probability
        #
        # |-------------|-------------------|-----------------|
        # |    K-mer    |   Ingroup = 15%   |   Other = 85%   |
        # |-------------|-------------------|-----------------|
        # |   PRESENT   |        80%        |        5%       |
        # |   ABSENT    |        20%        |       95%       |
        # |-------------|-------------------|-----------------|
        #
        # Ingroup = percentage os the sample records that belongs to ingroup.
        # Outgroup = percentage os the sample records that belongs to outgroup.
        #
        # ? --------------------------------------------------------------------

        sample_size = len(positive_reference) + len(negative_reference)
        kmers_post_probabilities: List[KmerBayesStats] = []
        available_kmers: Set = set()

        for reference in [*positive_reference, *negative_reference]:
            available_kmers.update(reference)

        for kmer in sorted(available_kmers):
            kmer_on_ingroup = __calculate_kmer_percentage_in_group(
                kmer,
                positive_reference,
            )

            kmer_on_outgroup = __calculate_kmer_percentage_in_group(
                kmer,
                negative_reference,
            )

            kmer_post_prob = __calculate_kmer_bayes_stats(
                __rescale_percent(len(positive_reference) * 100) / sample_size,
                __rescale_percent(kmer_on_ingroup.get_percentage()),
                __rescale_percent(kmer_on_outgroup.get_percentage()),
            )

            kmers_post_probabilities.append(
                KmerBayesStats(
                    kmer=kmer,
                    true_positives_post_prob=kmer_post_prob,
                )
            )

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(kmers_post_probabilities)

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))


def __rescale_percent(value: float) -> float:
    return value / 100


def __calculate_kmer_bayes_stats(
    ingroup_probability: float,
    ingroup_kmer_percentage: float,
    outgroup_kmer_percentage: float,
) -> float:
    # ? ------------------------------------------------------------------------
    # ? Validate args
    # ? ------------------------------------------------------------------------

    for arg in [
        ingroup_probability,
        ingroup_kmer_percentage,
        outgroup_kmer_percentage,
    ]:
        if arg > 1:
            raise

        if arg < 0:
            raise

    # ? ------------------------------------------------------------------------
    # ? Calculate the bayes stats
    #
    # |----------------|-------------------|-----------------|
    # |  K-mer = GAAT  |   Ingroup = 15%   |   Other = 85%   |
    # |----------------|-------------------|-----------------|
    # |   PRESENT      |        80%        |        5%       |
    # |   ABSENT       |        20%        |       95%       |
    # |----------------|-------------------|-----------------|
    #
    # Ingroup = percentage os the sample records that belongs to ingroup.
    # Outgroup = percentage os the sample records that belongs to outgroup.
    #
    # ---
    #
    # CALCULATION LOGIC:
    #
    # P(A|B) = Probability to belongs to Ingroup (A) given the k-mer existence
    # (B) being equals to 80%.
    #
    # P(B|A) = Probability of the k-mer exists (B) given the record belongs to
    # the Ingroup (A). These is a true positive (80%).
    #
    # P(A) = Probability to belongs to Ingroup (15%).
    #
    # P(AC) = Probability of record to NOT belongs to Ingroup (85%), in other
    # words, to belongs to the Outgroup.
    #
    # P(B|AC) = Probability to contains the k-mer (B) given the record not
    # belongs to the Ingroup (belonging to the Outgroup), being equals to 5%.
    #
    #                                  P(B|A) x P(A)
    #               P(A|B) =  ---------------------------------
    #                          P(B|A) x P(A) + P(B|AC) x P(AC)
    #
    #
    #                                   0.80 x 0.15
    #                 P(A|B) =  ---------------------------
    #                            0.80 x 0.15 + 0.05 x 0.85
    #
    #
    #                         P(A|B) = 0,74 = 74%
    # ? ------------------------------------------------------------------------

    true_positives = ingroup_kmer_percentage * ingroup_probability
    false_positives = outgroup_kmer_percentage * (1 - ingroup_probability)

    return true_positives / (true_positives + false_positives)


def __calculate_kmer_percentage_in_group(
    kmer: str,
    group_records: List[List[str]],
) -> __KmerStats:
    if len(group_records) < 5:
        raise

    kmer_on_group = 0

    for record in group_records:
        if kmer in record:
            kmer_on_group += 1

    return __KmerStats(
        group_length=len(group_records),
        records_on_group=kmer_on_group,
    )
