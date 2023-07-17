import clean_base.exceptions as c_exc
from attr import define, field
from clean_base.either import Either, left, right

from classeq.settings import LOGGER


@define(kw_only=True)
class __KmerStats:
    """Kmer stats.

    Attributes:
        group_length (int): The length of the group.
        records_on_group (int): The number of records on the group.

    """

    group_length: int = field()
    records_on_group: int = field()

    def get_percentage(self) -> float:
        return (self.records_on_group * 100) / self.group_length


@define(kw_only=True)
class KmerBayesStats:
    """Kmer Bayes stats.

    Attributes:
        kmer (str): The kmer.
        group_post_prob (float): The group posterior probability.

    """

    kmer: str
    group_post_prob: float


def calculate_joining_probability_for_group(
    prediction_target: list[str],
    positive_reference: list[list[str]],
    negative_reference: list[list[str]],
) -> Either[c_exc.MappedErrors, list[KmerBayesStats]]:
    """Calculate the probability of a sequence belongs to a clade.

    Args:
        prediction_target (list[str]): The sequence to be tested.
        positive_reference (list[list[str]]): The positive reference.
        negative_reference (list[list[str]]): The negative reference.

    Returns:
        Either[c_exc.MappedErrors, list[KmerBayesStats]]: The probability of
            the sequence belongs to the clade.

    Raises:
        c_exc.UseCaseError: If any of the arguments is not a list.
        c_exc.UseCaseError: If any of the sub-items of the arguments is
            not a list.

    """

    try:
        # ? --------------------------------------------------------------------
        # ? Validate args
        # ? --------------------------------------------------------------------

        if not isinstance(prediction_target, list):
            return left(
                c_exc.UseCaseError(
                    f"Prediction target is not a list: {prediction_target}"
                )
            )

        for name, group in [
            ("positive_reference", positive_reference),
            ("negative_reference", negative_reference),
        ]:
            if not isinstance(group, list):
                return left(
                    c_exc.UseCaseError(f"{name} is not a list: {group}")
                )

            for sub_item in group:
                if not isinstance(sub_item, list):
                    return left(
                        c_exc.UseCaseError(
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
        kmers_post_probabilities: list[KmerBayesStats] = []
        available_kmers: set[str] = set()

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
                    group_post_prob=kmer_post_prob,
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
    """Calculate the bayes stats.

    Args:
        ingroup_probability (float): The probability of the sequence belongs to
            the ingroup.
        ingroup_kmer_percentage (float): The percentage of the kmer on the
            ingroup.
        outgroup_kmer_percentage (float): The percentage of the kmer on the
            outgroup.

    Returns:
        float: The probability of the sequence belongs to the ingroup.

    Raises:
        c_exc.UseCaseError: If any of the arguments is not a float.
        c_exc.UseCaseError: If any of the arguments is not between 0
            and 1.

    """

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
    #
    # ? ------------------------------------------------------------------------

    true_positives = ingroup_kmer_percentage * ingroup_probability
    false_positives = outgroup_kmer_percentage * (1 - ingroup_probability)

    return true_positives / (true_positives + false_positives)


def __calculate_kmer_percentage_in_group(
    kmer: str,
    group_records: list[list[str]],
) -> __KmerStats:
    """Calculate the percentage of the k-mer on the group.

    Args:
        kmer (str): The k-mer to be calculated.
        group_records (list[list[str]]): The group of records to be
            calculated.

    Returns:
        __KmerStats: The k-mer stats.

    Raises:
        Exception: If the group has less than 5 records.

    """

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
