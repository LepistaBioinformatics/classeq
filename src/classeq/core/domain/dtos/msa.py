from pathlib import Path
from re import search
from typing import Any, DefaultDict, Self

from attr import define, field
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.msa_source_format import MsaSourceFormatEnum
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import (
    BASES,
    DEFAULT_KMER_SIZE,
    LOGGER,
    TEMP_INPUT_FILE_SUFFIX,
)


@define(kw_only=True)
class MsaSource:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    source_file_path: Path = field()
    file_format: MsaSourceFormatEnum = field()
    sequence_headers: list[int] = field()
    kmers_indices: KmersInverseIndices | None = field(default=None)

    # ? ------------------------------------------------------------------------
    # ? Public class methods
    # ? ------------------------------------------------------------------------

    @classmethod
    def from_dict(
        cls,
        content: dict[str, Any],
    ) -> Either[Self, c_exc.MappedErrors]:
        for key in [
            "source_file_path",
            "file_format",
            "sequence_headers",
            "kmers_indices",
        ]:
            if key not in content:
                return left(
                    c_exc.InvalidArgumentError(
                        f"Invalid content detected on parse `{MsaSource}`. "
                        f"{key}` key is empty.",
                        logger=LOGGER,
                    )
                )

        kmer_indices_either = KmersInverseIndices.from_dict(
            content=content.get("kmers_indices")
        )

        if kmer_indices_either.is_left:
            return kmer_indices_either

        return right(
            cls(
                source_file_path=Path(content.get("source_file_path")),  # type: ignore
                file_format=eval(content.get("file_format")),  # type: ignore
                sequence_headers=content.get("sequence_headers"),  # type: ignore
                kmers_indices=kmer_indices_either.value,
            )
        )

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

            sequence_headers: list[int] = list()

            cleaned_file_path = source_file_path.parent.joinpath(
                "".join(
                    [
                        source_file_path.stem,
                        ".",
                        TEMP_INPUT_FILE_SUFFIX,
                        source_file_path.suffix,
                    ]
                )
            )

            LOGGER.info(
                f"Sanitized MSA would be persisted to: {cleaned_file_path}"
            )

            with cleaned_file_path.open("w") as out:
                for record in SeqIO.parse(
                    handle=source_file_path,
                    format=format.value,
                ):
                    if cls.__check_sequence_sanity(record) is False:
                        LOGGER.warning(
                            f"Sequence `{record.id}` has invalid characters. "
                            + f"Only non redundant residuals ({', '.join(BASES)}) "
                            + "would be included in analysis."
                        )

                    sequence_headers.append(record.id)

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

    def initialize_kmer_indices(
        self,
        headers_map: DefaultDict[str, int],
    ) -> Either[bool, c_exc.MappedErrors]:
        try:
            indices_either = KmersInverseIndices.new(
                source_file_path=self.source_file_path,
                format=self.file_format,
                k_size=DEFAULT_KMER_SIZE,
                headers_map=headers_map,
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
