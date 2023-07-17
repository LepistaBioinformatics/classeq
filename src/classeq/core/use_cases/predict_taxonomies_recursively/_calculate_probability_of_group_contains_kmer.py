def calculate_probability_of_group_contains_kmer(
    prior: float,
    sequences_with_kmer: int,
    total_sequences: int,
) -> float:
    """Calculate the probability of a group contains a kmer.

    Description:
        This function calculates the probability of a group contains a kmer.
        The probability is calculated by the following formula:
        P(wi|G) = [ m(wi) + Pi ] / ( M + 1 ) where m(wi) is the number of the
        group sequences containing word wi, Pi is the prior probability for the
        kmer in the overall training dataset, and M is the dataset size for the
        target training group.

    Args:
        prior (float): The prior probability for the kmer in the overall
            training dataset.
        sequences_with_kmer (int): The number of the group sequences
            containing word wi.
        total_sequences (int): The dataset size for the target training group.

    Returns:
        float: The probability of a group contains a kmer.

    Raises:
        Exception: If `total_sequences` is greater than `sequences_with_kmer`.

    """

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
