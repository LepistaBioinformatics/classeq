from collections import defaultdict
from pathlib import Path
from re import search
from typing import DefaultDict, List, Self

from attr import define, field
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.kmer_inverse_index import KmerInverseIndex
from classeq.core.domain.dtos.msa_source_format import MsaSourceFormatEnum
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import BASES, DEFAULT_KMER_SIZE, LOGGER


@define(kw_only=True)
class MsaSource:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    source_file_path: Path = field()
    file_format: MsaSourceFormatEnum = field()
    sequence_headers: DefaultDict[str, int] = field()
    raw_sequences: List[SeqRecord] | None = field(init=False, default=None)
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

            for index, record in enumerate(
                SeqIO.parse(
                    handle=source_file_path,
                    format=format.value,
                )
            ):
                if cls.__check_sequence_sanity(record) is False:
                    return left(
                        c_exc.InvalidArgumentError(
                            f"Invalid sequence: {record.id}"
                        )
                    )

                sequence_headers[record.id] = index

            return right(
                cls(
                    source_file_path=source_file_path,
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
                        f"Unexpected error on generate kmer indices.",
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
    def __check_sequence_sanity(sequence: Seq) -> bool:
        return search(f"[^{''.join(BASES)}]", str(sequence.seq).upper()) is True
