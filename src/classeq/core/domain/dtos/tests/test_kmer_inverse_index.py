import logging
from collections import defaultdict
from os import getenv
from pathlib import Path
from typing import DefaultDict
from unittest import TestCase

from Bio import SeqIO

from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.msa import MsaSourceFormatEnum


class KmerInverseIndexTest(TestCase):
    def setUp(self) -> None:
        logging.basicConfig(level=logging.DEBUG)
        self.__logger = logging.getLogger()

    def test_if_it_work(self) -> None:
        source_file_path = Path(str(getenv("MOCK_MSA")))
        format = MsaSourceFormatEnum.FASTA

        sequence_headers: DefaultDict[str, int] = defaultdict()

        for index, record in enumerate(
            SeqIO.parse(
                handle=source_file_path,
                format=format.value,
            )
        ):
            sequence_headers[record.id] = index

        indices = KmersInverseIndices.new(
            source_file_path=source_file_path,
            format=MsaSourceFormatEnum.FASTA,
            k_size=8,
            headers_map=sequence_headers,
        )

        self.assertFalse(indices.is_left)
        self.assertTrue(indices.is_right)
        self.assertIsInstance(indices.value, KmersInverseIndices)


if __name__ == "__main__":
    from unittest import main

    main()
