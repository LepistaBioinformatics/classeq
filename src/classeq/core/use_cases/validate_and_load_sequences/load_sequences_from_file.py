from pathlib import Path

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.msa import MsaSource, MsaSourceFormatEnum
from classeq.core.domain.utils.either import Either, left
from classeq.settings import LOGGER


def load_sequences_from_file(
    source_file_path: Path,
    format: MsaSourceFormatEnum,
) -> Either[MsaSource, c_exc.MappedErrors]:
    try:
        return MsaSource.new(
            source_file_path=source_file_path,
            format=format,
        )

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))
