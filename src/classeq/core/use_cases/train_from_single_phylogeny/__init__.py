import gzip
from json import dump

from attrs import asdict

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.priors import TreePriors
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER

from .estimate_clade_specific_priors import estimate_clade_specific_priors


def train_from_single_phylogeny(
    references: ReferenceSet,
) -> Either[bool, c_exc.MappedErrors]:
    try:
        # ? --------------------------------------------------------------------
        # ? Validate args
        # ? --------------------------------------------------------------------

        if not isinstance(references, ReferenceSet):
            return left(
                c_exc.InvalidArgumentError(
                    f"Argument `references` should be a `{ReferenceSet}` instance.",
                    exp=True,
                    logger=LOGGER,
                )
            )

        # ? --------------------------------------------------------------------
        # ? Recursively train for clades
        # ? --------------------------------------------------------------------

        train_response_either = estimate_clade_specific_priors(
            references=references
        )

        if train_response_either.is_left:
            return train_response_either

        train_response: TreePriors = train_response_either.value

        # ? --------------------------------------------------------------------
        # ? Persist results to file
        # ? --------------------------------------------------------------------

        tree_source = references.tree.source_file_path
        train_output_file_path = tree_source.parent.joinpath(
            ".".join(
                [
                    tree_source.stem,
                    "train-source",
                    "json",
                    "gz",
                ]
            )
        )

        LOGGER.info(
            f"Train output file would be persisted to: {train_output_file_path}"
        )

        with gzip.open(
            train_output_file_path, "wt", encoding="utf-8"
        ) as out_gz:
            dump(
                asdict(train_response),
                out_gz,
                indent=4,
                default=str,
                sort_keys=True,
            )

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(True)

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))
