from pathlib import Path

import clean_base.exceptions as c_exc
from clean_base.either import Either, left

from classeq.core.domain.dtos.msa import MsaSource, MsaSourceFormatEnum
from classeq.settings import LOGGER


def load_and_sanitize_sequences(
    source_file_path: Path,
    format: MsaSourceFormatEnum,
) -> Either[c_exc.MappedErrors, MsaSource]:
    try:
        return MsaSource.new(
            source_file_path=source_file_path,
            format=format,
        )

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))
