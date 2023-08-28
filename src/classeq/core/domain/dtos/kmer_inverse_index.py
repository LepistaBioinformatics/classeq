from collections import defaultdict
from hashlib import md5
from itertools import islice
from pathlib import Path
from re import search
from typing import Any, Generator, Iterator, Self

import clean_base.exceptions as c_exc
from attrs import define, field, frozen
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from clean_base.either import Either, right

from classeq.core.domain.dtos.msa_source_format import MsaSourceFormatEnum
from classeq.core.domain.dtos.ordered_tuple import OrderedTuple
from classeq.core.domain.dtos.strand import StrandEnum
from classeq.settings import BASES, LOGGER


@define(kw_only=True)
class KmerIndex:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    kmer: str = field()
    records: OrderedTuple = field()

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

    def to_dict(self) -> dict[str, Any]:
        return {
            "kmer": self.kmer,
            "records": self.records,
        }

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

        return right(cls(kmer=kmer, records=OrderedTuple(records)))


@frozen(kw_only=True)
class KmersInverseIndices:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    indices: tuple[KmerIndex, ...] = field(default=tuple())

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return {
            "indices": [i.to_dict() for i in self.indices],
        }

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

        return right(cls(indices=tuple(kmer_indices)))

    @classmethod
    def new(
        cls,
        source_file_path: Path,
        format: MsaSourceFormatEnum,
        k_size: int,
        headers_map: defaultdict[str, int],
        strand: StrandEnum,
    ) -> Either[c_exc.MappedErrors, Self]:
        try:
            if not source_file_path.is_file():
                return c_exc.DadaTransferObjectError(
                    f"Invalid path: {source_file_path}"
                )()

            kmer_indices: defaultdict[str, set[int]] = defaultdict(set)

            for record in SeqIO.parse(
                handle=source_file_path,
                format=format.value,
            ):
                for kmer in cls.generate_kmers(
                    dna_sequence=cls.sanitize_sequence(record, True),
                    k_size=k_size,
                    strand=strand,
                ):
                    if (header_index := headers_map.get(record.id)) is None:
                        return c_exc.DadaTransferObjectError(
                            "Unexpected unmatch between kmer indices and "
                            + f"MSA sequence headers: {record.id}",
                            exp=True,
                            logger=LOGGER,
                        )()

                    kmer_indices[kmer].add(header_index)

            return right(
                cls(
                    indices=tuple(
                        sorted(
                            [
                                KmerIndex(
                                    kmer=kmer,
                                    records=OrderedTuple(codes),
                                )
                                for kmer, codes in kmer_indices.items()
                            ],
                            key=lambda i: i.__hash__(),
                        )
                    )
                )
            )

        except Exception as exc:
            return c_exc.CreationError(exc, logger=LOGGER)()

    @classmethod
    def generate_kmers(
        cls,
        dna_sequence: Seq,
        k_size: int,
        strand: StrandEnum,
    ) -> Generator[str, None, None]:
        if strand in [StrandEnum.PLUS, StrandEnum.BOTH]:
            for kmer in cls.__generate_kmers(
                dna_sequence=dna_sequence,
                k_size=k_size,
            ):
                yield kmer

        if strand in [StrandEnum.MINUS, StrandEnum.BOTH]:
            for kmer in cls.__generate_kmers(
                dna_sequence=dna_sequence.reverse_complement(),
                k_size=k_size,
            ):
                yield kmer

    # ? ------------------------------------------------------------------------
    # ? Public static methods
    # ? ------------------------------------------------------------------------

    @staticmethod
    def check_sequence_sanity(sequence: SeqRecord) -> bool:
        return search(f"[^{''.join(BASES)}]", str(sequence.seq).upper()) is None

    @staticmethod
    def sanitize_sequence(
        sequence: SeqRecord,
        as_seq: bool = False,
    ) -> str:
        seq = "".join(
            [
                letter
                for letter in Seq(
                    sequence.seq.replace("-", "").replace("?", "")
                ).upper()
                if letter in BASES
            ]
        )

        if as_seq is True:
            return Seq(seq)
        return seq

    # ? ------------------------------------------------------------------------
    # ? Private static methods
    # ? ------------------------------------------------------------------------

    @staticmethod
    def __generate_kmers(
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
