from collections import defaultdict
from pathlib import Path
from typing import Any, Self

import clean_base.exceptions as c_exc
from attr import define, field
from clean_base.either import Either, right

from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.msa_source_format import MsaSourceFormatEnum
from classeq.core.domain.dtos.strand import StrandEnum
from classeq.settings import LOGGER


@define(kw_only=True)
class MsaSource:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

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

        if (kmers_indices := content.get("kmers_indices")) is None:
            return c_exc.DadaTransferObjectError(
                "Invalid content detected on parse `kmers_indices` from JSON dump",
                logger=LOGGER,
            )()

        kmer_indices_either = KmersInverseIndices.from_dict(
            content=kmers_indices
        )

        if kmer_indices_either.is_left:
            return kmer_indices_either

        return right(
            cls(
                file_format=eval(content.get("file_format")),  # type: ignore
                sequence_headers=content.get("sequence_headers"),  # type: ignore
                kmers_indices=kmer_indices_either.value,
            )
        )

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return {
            "file_format": self.file_format,
            "sequence_headers": self.sequence_headers,
            "kmers_indices": (
                self.kmers_indices.to_dict() if self.kmers_indices else None
            ),
        }

    def initialize_kmer_indices(
        self,
        source_file_path: Path,
        headers_map: defaultdict[str, int],
        k_size: int,
        strand: StrandEnum,
    ) -> Either[c_exc.MappedErrors, bool]:
        try:
            if (
                indices_either := KmersInverseIndices.new(
                    source_file_path=source_file_path,
                    format=self.file_format,
                    k_size=k_size,
                    headers_map=headers_map,
                    strand=strand,
                )
            ).is_left:
                return c_exc.DadaTransferObjectError(
                    "Unexpected error on generate kmer indices.",
                    prev=indices_either.value,
                    logger=LOGGER,
                )()

            self.kmers_indices = indices_either.value

            return right(True)

        except Exception as exc:
            return c_exc.CreationError(exc, logger=LOGGER)()
