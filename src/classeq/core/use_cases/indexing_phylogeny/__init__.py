import tarfile
from tempfile import TemporaryDirectory
import gzip
from json import dump
from pathlib import Path

import clean_base.exceptions as c_exc
from attrs import asdict
from clean_base.either import Either, right

from classeq.core.domain.dtos.priors import TreePriors
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.settings import (
    LOGGER,
    MINIMUM_CLADE_SIZE,
    TRAIN_SOURCE_OUTPUT_FILE_NAME,
)

from ._fetch_clade_specific_kmers import fetch_clade_specific_kmers


def indexing_phylogeny(
    references: ReferenceSet,
    indexing_file_path: Path,
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

        LOGGER.info("Output file:")
        LOGGER.info(
            f"\t{indexing_file_path.relative_to(indexing_file_path.parent)}"
        )

        tmp_dir = TemporaryDirectory()

        train_output_file_path = Path(tmp_dir.name).joinpath(
            TRAIN_SOURCE_OUTPUT_FILE_NAME
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

        with tarfile.open(indexing_file_path, "a:") as tar:
            tar.add(
                train_output_file_path,
                arcname=train_output_file_path.name,
            )

        tmp_dir.cleanup()

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(True)

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
