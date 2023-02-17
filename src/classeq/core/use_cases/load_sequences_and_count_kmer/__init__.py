from pathlib import Path
from classeq.core.domain.dtos.msa import MsaSourceFormatEnum

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER

from .load_sequences_from_file import load_sequences_from_file


def load_sequences_and_count_kmer(
    source_file_path: Path,
    format: MsaSourceFormatEnum,
) -> Either[bool, c_exc.MappedErrors]:
    try:
        # ? --------------------------------------------------------------------
        # ? Validate args
        # ? --------------------------------------------------------------------

        if not source_file_path.is_file():
            return left(
                c_exc.InvalidArgumentError(
                    f"Input file not exists: {source_file_path}"
                )
            )

        # ? --------------------------------------------------------------------
        # ? Load input file content
        # ? --------------------------------------------------------------------

        sequences_either: Either = load_sequences_from_file(
            source_file_path=source_file_path,
            format=format,
        )

        if sequences_either.is_left:
            return left(
                c_exc.UseCaseError(
                    "Unexpected error on load sequences.",
                    prev=sequences_either.value,
                    logger=LOGGER,
                )
            )

        sequences = sequences_either.value

        LOGGER.debug(f"sequences: {sequences}")

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(True)

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))
