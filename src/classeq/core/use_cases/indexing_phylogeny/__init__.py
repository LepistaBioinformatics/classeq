import gzip
from json import dump

import clean_base.exceptions as c_exc
from attrs import asdict
from clean_base.either import Either, right

from classeq.core.domain.dtos.priors import TreePriors
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.settings import LOGGER, MINIMUM_CLADE_SIZE

from ._fetch_clade_specific_kmers import fetch_clade_specific_kmers


def indexing_phylogeny(
    references: ReferenceSet,
    min_clade_size: int = MINIMUM_CLADE_SIZE,
) -> Either[c_exc.MappedErrors, bool]:
    try:
        # ? --------------------------------------------------------------------
        # ? Validate args
        # ? --------------------------------------------------------------------

        if not isinstance(references, ReferenceSet):
            return c_exc.UseCaseError(
                f"Argument `references` should be a `{ReferenceSet}` instance.",
                exp=True,
                logger=LOGGER,
            )()

        # ? --------------------------------------------------------------------
        # ? Recursively train for clades
        # ? --------------------------------------------------------------------

        if (
            train_response_either := fetch_clade_specific_kmers(
                references=references,
                min_clade_size=min_clade_size,
            )
        ).is_left:
            return train_response_either

        train_response: TreePriors = train_response_either.value

        # ? --------------------------------------------------------------------
        # ? Persist results to file
        # ? --------------------------------------------------------------------

        tree_source = references.tree.newick_file_path
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

        LOGGER.info("Train output file would be persisted to:")
        LOGGER.info(f"\t{train_output_file_path}")

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
        return c_exc.UseCaseError(exc, logger=LOGGER)()
