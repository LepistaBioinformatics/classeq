from hashlib import md5
from collections import defaultdict
from itertools import islice
from pathlib import Path
from typing import DefaultDict, Iterator, Self, Set  # type: ignore

from attrs import define, field, frozen
from Bio import SeqIO

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.msa_source_format import MsaSourceFormatEnum
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER


@define(kw_only=True)
class KmerIndex:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    kmer: str = field()
    records: Set[int] = field(default=set())

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
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def contains(
        self,
        target: int,
    ) -> bool:
        records = sorted(self.records)
        first = 0
        last = len(records) - 1
        index: int | None = None

        while (first <= last) and (index is None):
            mid = (first + last) // 2

            if records[mid] == target:
                index = mid
            else:
                if target < records[mid]:
                    last = mid - 1
                else:
                    first = mid + 1

        return index is None


@frozen(kw_only=True)
class KmerInverseIndex:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    indices: Set[KmerIndex] = field(default=set())
    hashes: Set[int] = field()

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
    ) -> Either[Self, c_exc.MappedErrors]:
        try:
            if not source_file_path.is_file():
                return left(
                    c_exc.InvalidArgumentError(
                        f"Invalid path: {source_file_path}"
                    )
                )

            kmer_indices: DefaultDict[str, Set[int]] = defaultdict(set)

            for record in SeqIO.parse(
                handle=source_file_path,
                format=format.value,
            ):
                for kmer in set(
                    [
                        kmer
                        for kmer in cls.__kmer_gen(
                            dna_sequence=record.seq,
                            k_size=k_size,
                        )
                    ]
                ):
                    header_index = headers_map.get(record.id)

                    if header_index is None:
                        return left(
                            c_exc.CreationError(
                                "Unexpected unmatch between kmer indices and "
                                + f"MSA sequence headers: {record.id}",
                                exp=True,
                                logger=LOGGER,
                            )
                        )

                    kmer_indices[kmer].add(header_index)

            indices = set(
                KmerIndex(
                    kmer=kmer,
                    records=indices,
                )
                for kmer, indices in kmer_indices.items()
            )

            return right(
                cls(
                    indices=indices,
                    hashes=set(i.__hash__() for i in indices),
                )
            )

        except Exception as exc:
            return left(c_exc.CreationError(exc, logger=LOGGER))

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def index_of(
        self,
        kmer: str,
    ) -> int | None:
        hashed_kmer = int(md5(kmer.encode("utf-8")).hexdigest(), 16)
        records = sorted(self.hashes)
        first = 0
        last = len(records) - 1
        index: int | None = None

        while (first <= last) and (index is None):
            mid = (first + last) // 2

            if records[mid] == hashed_kmer:
                index = mid
            else:
                if hashed_kmer < records[mid]:
                    last = mid - 1
                else:
                    first = mid + 1

        return index

    # ? ------------------------------------------------------------------------
    # ? Private static methods
    # ? ------------------------------------------------------------------------

    @staticmethod
    def __kmer_gen(
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
