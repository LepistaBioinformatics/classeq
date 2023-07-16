from __future__ import annotations

from collections import defaultdict
from hashlib import md5
from itertools import islice
from pathlib import Path
from typing import Any, DefaultDict, Iterator, Self

import clean_base.exceptions as c_exc
from attrs import define, field, frozen
from Bio import SeqIO
from clean_base.either import Either, right

from classeq.core.domain.dtos.msa_source_format import MsaSourceFormatEnum
from classeq.settings import LOGGER


class OrderedKmerRecords(tuple[int, ...]):
    def __new__(cls, values: tuple[int, ...] | list[int] | set[int]) -> Self:
        return super(OrderedKmerRecords, cls).__new__(cls, tuple(sorted(values)))  # type: ignore


@define(kw_only=True)
class KmerIndex:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    kmer: str = field()
    records: OrderedKmerRecords = field()

    # ? ------------------------------------------------------------------------
    # ? Life cycle hook methods
    # ? ------------------------------------------------------------------------

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, KmerIndex):
            return NotImplemented
        return self.__hash__() == other.__hash__()

    def __ne__(self, other: object) -> bool:
        # Not strictly necessary, but to avoid having both x==y and x!=y True at
        # the same time.
        return not (self.__hash__() == other.__hash__())

    def __hash__(self) -> int:
        return int(md5(self.kmer.encode("utf-8")).hexdigest(), 16)

    # ? ------------------------------------------------------------------------
    # ? Public class methods
    # ? ------------------------------------------------------------------------

    @classmethod
    def from_dict(
        cls,
        content: dict[str, Any],
    ) -> Either[c_exc.MappedErrors, Self]:
        if (kmer := content.get("kmer")) is None:
            return c_exc.DadaTransferObjectError(
                "Invalid content detected on parse `kmer` from JSON dump",
                logger=LOGGER,
            )()

        if (records := content.get("records")) is None:
            return c_exc.DadaTransferObjectError(
                "Invalid content detected on load `records` from JSON dump",
                logger=LOGGER,
            )()

        return right(cls(kmer=kmer, records=OrderedKmerRecords(records)))

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def contains(
        self,
        target: int,
    ) -> int | None:
        first = 0
        last = len(self.records) - 1
        index: int | None = None

        while (first <= last) and (index is None):
            mid = (first + last) // 2

            if self.records[mid] == target:
                return mid
            else:
                if target < self.records[mid]:
                    last = mid - 1
                else:
                    first = mid + 1

        return index


@frozen(kw_only=True)
class KmersInverseIndices:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    indices: tuple[KmerIndex, ...] = field(default=tuple())
    hashes: tuple[int, ...] = field()

    # ? ------------------------------------------------------------------------
    # ? Public class methods
    # ? ------------------------------------------------------------------------

    @classmethod
    def from_dict(
        cls,
        content: dict[str, Any],
    ) -> Either[c_exc.MappedErrors, Self]:
        for key in [
            "indices",
            "hashes",
        ]:
            if key not in content:
                return c_exc.DadaTransferObjectError(
                    f"Invalid content detected on parse `{KmersInverseIndices}`. "
                    f"{key}` key is empty.",
                    logger=LOGGER,
                )()

        kmer_indices: list[KmerIndex] = []
        for index in content.get("indices"):  # type: ignore
            kmer_index_either = KmerIndex.from_dict(content=index)

            if kmer_index_either.is_left:
                return kmer_index_either

            kmer_indices.append(kmer_index_either.value)

        return right(
            cls(
                indices=tuple(kmer_indices),
                hashes=tuple(sorted(content.get("hashes"))),  # type: ignore
            )
        )

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    @classmethod
    def new(
        cls,
        source_file_path: Path,
        format: MsaSourceFormatEnum,
        k_size: int,
        headers_map: DefaultDict[str, int],
    ) -> Either[c_exc.MappedErrors, Self]:
        try:
            if not source_file_path.is_file():
                return c_exc.DadaTransferObjectError(
                    f"Invalid path: {source_file_path}"
                )()

            kmer_indices: DefaultDict[str, set[int]] = defaultdict(set)

            for record in SeqIO.parse(
                handle=source_file_path,
                format=format.value,
            ):
                for kmer in set(
                    [
                        kmer
                        for kmer in cls.generate_kmers(
                            dna_sequence=record.seq,
                            k_size=k_size,
                        )
                    ]
                ):
                    header_index = headers_map.get(record.id)

                    if header_index is None:
                        return c_exc.DadaTransferObjectError(
                            "Unexpected unmatch between kmer indices and "
                            + f"MSA sequence headers: {record.id}",
                            exp=True,
                            logger=LOGGER,
                        )()

                    kmer_indices[kmer].add(header_index)

            indices = tuple(
                sorted(
                    [
                        KmerIndex(
                            kmer=kmer,
                            records=OrderedKmerRecords(codes),
                        )
                        for kmer, codes in kmer_indices.items()
                    ],
                    key=lambda i: i.__hash__(),
                )
            )

            return right(
                cls(
                    indices=indices,
                    hashes=tuple(i.__hash__() for i in indices),
                )
            )

        except Exception as exc:
            return c_exc.CreationError(exc, logger=LOGGER)()

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def index_of(
        self,
        kmer: str,
    ) -> int | None:
        # ? Initialize search params
        first = 0
        last = len(self.hashes) - 1
        index: int | None = None

        # ? Convert the kmer to a md5 representation
        hashed_kmer = int(md5(kmer.encode("utf-8")).hexdigest(), 16)

        while (first <= last) and (index is None):
            mid = (first + last) // 2

            if self.hashes[mid] == hashed_kmer:
                return mid
            else:
                if hashed_kmer < self.hashes[mid]:
                    last = mid - 1
                else:
                    first = mid + 1

        return index

    # ? ------------------------------------------------------------------------
    # ? Public static methods
    # ? ------------------------------------------------------------------------

    @staticmethod
    def generate_kmers(
        dna_sequence: str,
        k_size: int,
    ) -> Iterator[str]:
        iter_sequence = iter(dna_sequence)
        iter_content = tuple(islice(iter_sequence, k_size))

        if len(iter_content) == k_size:
            yield "".join(iter_content)

        for ele in iter_sequence:
            iter_content = iter_content[1:] + (ele,)
            yield "".join(iter_content)
