from __future__ import annotations

from enum import Enum
from pathlib import Path
from re import search
from typing import List, Optional

from attr import define, field
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER


class MsaSourceFormatEnum(Enum):
    FASTA = "fasta"


@define(kw_only=True)
class MsaSource:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    source_file_path: Path = field()
    fasta_headers: List[str] = field()
    raw_sequences: Optional[List[SeqRecord]] = field()

    # ? ------------------------------------------------------------------------
    # ? Public class methods
    # ? ------------------------------------------------------------------------

    @classmethod
    def new(
        cls,
        source_file_path: Path,
        format: MsaSourceFormatEnum,
    ) -> Either[MsaSource, c_exc.MappedErrors]:
        try:
            if not source_file_path.is_file():
                return left(
                    c_exc.InvalidArgumentError(
                        f"Invalid path: {source_file_path}"
                    )
                )

            fasta_headers: List[str] = []

            for record in SeqIO.parse(
                handle=source_file_path,
                format=format.value,
            ):
                if cls.__check_sequence_sanity(record) is False:
                    return left(
                        c_exc.InvalidArgumentError(
                            f"Invalid sequence: {record.id}"
                        )
                    )

                fasta_headers.append(record.id)

            return right(
                cls(  # type: ignore
                    source_file_path=source_file_path,
                    fasta_headers=fasta_headers,
                    raw_sequences=None,
                )
            )

        except Exception as exc:
            return left(c_exc.UseCaseError(exc, logger=LOGGER))

    # ? ------------------------------------------------------------------------
    # ? Private instance methods
    # ? ------------------------------------------------------------------------

    @staticmethod
    def __check_sequence_sanity(sequence: Seq) -> bool:
        return search("[^ATCG]", str(sequence.seq).upper()) is True
