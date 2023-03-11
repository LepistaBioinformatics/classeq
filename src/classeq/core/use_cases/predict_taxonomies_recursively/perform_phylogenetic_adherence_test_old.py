import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER

from .calculate_joining_probability_for_group import (
    calculate_joining_probability_for_group,
)


# !
# ! Usar os resultados de input, sister e outgroup para estimar com bayesiana:
# ! https://towardsdatascience.com/estimating-probabilities-with-bayesian-modeling-in-python-7144be007815
# !


def perform_phylogenetic_adherence_test(
    prediction_target: list[str],
    ingroup: list[list[str]],
    sister: list[list[str]],
    outgroup: list[list[str]],
) -> Either[c_exc.MappedErrors, bool]:
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
            ("ingroup", ingroup),
            ("sister", sister),
            ("outgroup", outgroup),
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
        # ? Generate priors by group
        # ? --------------------------------------------------------------------

        ingroup_pp = calculate_joining_probability_for_group(
            prediction_target=prediction_target,
            positive_reference=ingroup,
            negative_reference=[*sister, *outgroup],
        )

        for item in ingroup_pp.value:
            LOGGER.debug(f"ingroup_pp: {item}")

        sister_pp = calculate_joining_probability_for_group(
            prediction_target=prediction_target,
            positive_reference=sister,
            negative_reference=[*outgroup, *ingroup],
        )

        for item in sister_pp.value:
            LOGGER.debug(f"sister_pp: {item}")

        outgroup_pp = calculate_joining_probability_for_group(
            prediction_target=prediction_target,
            positive_reference=outgroup,
            negative_reference=[*ingroup, *sister],
        )

        for item in outgroup_pp.value:
            LOGGER.debug(f"outgroup_pp: {item}")

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(True)

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))
