from pathlib import Path
from re import search
from typing import Any, DefaultDict, Self

import clean_base.exceptions as c_exc
from attr import define, field
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from clean_base.either import Either, right

from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.msa_source_format import MsaSourceFormatEnum
from classeq.settings import (
    BASES,
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
    ) -> Either[c_exc.MappedErrors, Self]:
        for key in [
            "source_file_path",
            "file_format",
            "sequence_headers",
            "kmers_indices",
        ]:
            if key not in content:
                return c_exc.DadaTransferObjectError(
                    f"Invalid content detected on parse `{MsaSource}`. "
                    f"{key}` key is empty.",
                    logger=LOGGER,
                )()

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
        output_directory: Path,
    ) -> Either[c_exc.MappedErrors, Self]:
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

            LOGGER.info("Sanitized MSA should be persisted to:")
            LOGGER.info(f"\t{cleaned_file_path.relative_to(Path.cwd())}")

            with cleaned_file_path.open("w") as out:
                for record in SeqIO.parse(
                    handle=source_file_path,
                    format=format.value,
                ):
                    if cls.__check_sequence_sanity(record) is False:
                        LOGGER.warning(
                            f"Sequence `{record.id}` has invalid characters. "
                            + f"Only non redundant residuals ({', '.join(BASES)}) "
                            + "should be included in analysis."
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
            return c_exc.CreationError(exc, logger=LOGGER)()

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def initialize_kmer_indices(
        self,
        headers_map: DefaultDict[str, int],
        k_size: int,
    ) -> Either[c_exc.MappedErrors, bool]:
        try:
            indices_either = KmersInverseIndices.new(
                source_file_path=self.source_file_path,
                format=self.file_format,
                k_size=k_size,
                headers_map=headers_map,
            )

            if indices_either.is_left:
                return c_exc.DadaTransferObjectError(
                    "Unexpected error on generate kmer indices.",
                    prev=indices_either.value,
                    logger=LOGGER,
                )()

            self.kmers_indices = indices_either.value

            return right(True)

        except Exception as exc:
            return c_exc.CreationError(exc, logger=LOGGER)()

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
