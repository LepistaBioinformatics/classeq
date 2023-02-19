from collections import defaultdict
from pathlib import Path
from re import search
from typing import DefaultDict, Self  # type: ignore

from attr import define, field
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.kmer_inverse_index import KmerInverseIndex
from classeq.core.domain.dtos.msa_source_format import MsaSourceFormatEnum
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import (
    BASES,
    DEFAULT_KMER_SIZE,
    LOGGER,
    TEMP_INPUT_FILE_PREFIX,
)


@define(kw_only=True)
class MsaSource:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    source_file_path: Path = field()
    file_format: MsaSourceFormatEnum = field()
    sequence_headers: DefaultDict[str, int] = field()
    kmers_indices: KmerInverseIndex | None = field(init=False, default=None)

    # ? ------------------------------------------------------------------------
    # ? Public class methods
    # ? ------------------------------------------------------------------------

    @classmethod
    def new(
        cls,
        source_file_path: Path,
        format: MsaSourceFormatEnum,
    ) -> Either[Self, c_exc.MappedErrors]:
        try:
            if not source_file_path.is_file():
                return left(
                    c_exc.InvalidArgumentError(
                        f"Invalid path: {source_file_path}"
                    )
                )

            sequence_headers: DefaultDict[str, int] = defaultdict()

            cleaned_file_path = source_file_path.parent.joinpath(
                "".join(
                    [
                        source_file_path.stem,
                        ".",
                        TEMP_INPUT_FILE_PREFIX,
                        source_file_path.suffix,
                    ]
                )
            )

            LOGGER.info(
                f"Sanitized file would be persisted to: {cleaned_file_path}"
            )

            with cleaned_file_path.open("w") as out:
                for index, record in enumerate(
                    SeqIO.parse(
                        handle=source_file_path,
                        format=format.value,
                    )
                ):
                    if cls.__check_sequence_sanity(record) is False:
                        LOGGER.warning(
                            f"Sequence `{record.id}` has invalid characters. "
                            + f"Only non redundant residuals ({', '.join(BASES)}) "
                            + "would be included in analysis."
                        )

                    sequence_headers[record.id] = index

                    sanitized_sequence = cls.__sanitize_sequence(record)
                    out.write(f">{record.id}\n{sanitized_sequence}\n")

            return right(
                cls(
                    source_file_path=cleaned_file_path,
                    file_format=format,
                    sequence_headers=sequence_headers,
                )
            )

        except Exception as exc:
            return left(c_exc.CreationError(exc, logger=LOGGER))

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def initialize_kmer_indices(self) -> Either[bool, c_exc.MappedErrors]:
        try:
            indices_either = KmerInverseIndex.new(
                self.source_file_path,
                self.file_format,
                DEFAULT_KMER_SIZE,
                self.sequence_headers,
            )

            if indices_either.is_left:
                return left(
                    c_exc.InvalidArgumentError(
                        "Unexpected error on generate kmer indices.",
                        prev=indices_either.value,
                        logger=LOGGER,
                    )
                )

            self.kmers_indices = indices_either.value

            return right(True)

        except Exception as exc:
            return left(c_exc.CreationError(exc, logger=LOGGER))

    # ? ------------------------------------------------------------------------
    # ? Private instance methods
    # ? ------------------------------------------------------------------------

    @staticmethod
    def __check_sequence_sanity(sequence: SeqRecord) -> bool:
        return search(f"[^{''.join(BASES)}]", str(sequence.seq).upper()) is None

    @staticmethod
    def __sanitize_sequence(sequence: SeqRecord) -> str:
        return "".join(
            [letter for letter in str(sequence.seq).upper() if letter in BASES]
        )