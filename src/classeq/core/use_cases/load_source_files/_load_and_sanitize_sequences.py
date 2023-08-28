from pathlib import Path

import clean_base.exceptions as c_exc
from Bio import SeqIO
from clean_base.either import Either, left, right
from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices

from classeq.core.domain.dtos.msa import MsaSource, MsaSourceFormatEnum
from classeq.settings import BASES, LOGGER, TEMP_INPUT_FILE_SUFFIX


def load_and_sanitize_sequences(
    source_file_path: Path,
    format: MsaSourceFormatEnum,
    output_directory: Path,
) -> Either[c_exc.MappedErrors, tuple[Path, MsaSource]]:
    try:
        if not source_file_path.is_file():
            return c_exc.DadaTransferObjectError(
                f"Invalid path: {source_file_path}"
            )()

        sequence_headers: list[int] = list()

        cleaned_file_path = output_directory.joinpath(
            "".join(
                [
                    source_file_path.stem,
                    ".",
                    TEMP_INPUT_FILE_SUFFIX,
                    source_file_path.suffix,
                ]
            )
        )

        LOGGER.info("Sanitized MSA file:")
        LOGGER.info(f"\t{cleaned_file_path.relative_to(output_directory)}")

        with cleaned_file_path.open("w+") as out:
            for record in SeqIO.parse(
                handle=source_file_path,
                format=format.value,
            ):
                record.seq = KmersInverseIndices.sanitize_sequence(record, True)
                if KmersInverseIndices.check_sequence_sanity(record) is False:
                    LOGGER.warning(
                        f"Sequence `{record.id}` has invalid characters. "
                        + f"Only non redundant residuals ({', '.join(BASES)}) "
                        + "should be included in analysis."
                    )

                sequence_headers.append(record.id)
                out.write(f">{record.id}\n{record.seq}\n")

        return right(
            (
                cleaned_file_path,
                MsaSource(
                    file_format=format,
                    sequence_headers=sequence_headers,
                ),
            )
        )

    except Exception as exc:
        return left(c_exc.UseCaseError(exc, logger=LOGGER))
