import logging
from collections import defaultdict
from os import getenv
from pathlib import Path
from typing import DefaultDict
from unittest import TestCase

from classeq.core.domain.dtos.msa import MsaSource
from classeq.core.domain.dtos.msa_source_format import MsaSourceFormatEnum


class MsaTest(TestCase):
    def setUp(self) -> None:
        logging.basicConfig(level=logging.DEBUG)
        self.__logger = logging.getLogger()

    def test_if_it_work(self) -> None:
        source_file_path = Path(str(getenv("MOCK_MSA")))
        format = MsaSourceFormatEnum.FASTA

        msa_source_either = MsaSource.new(
            source_file_path=source_file_path,
            format=format,
        )

        self.assertFalse(msa_source_either.is_left)
        self.assertTrue(msa_source_either.is_right)
        self.assertIsInstance(msa_source_either.value, MsaSource)

        msa_source: MsaSource = msa_source_either.value

        labels_map: DefaultDict[str, int] = defaultdict()

        for index, header in enumerate(msa_source.sequence_headers):
            labels_map[header] = index

        bool_response_either = msa_source.initialize_kmer_indices(
            headers_map=labels_map
        )

        self.assertFalse(bool_response_either.is_left)
        self.assertTrue(bool_response_either.is_right)
        self.assertTrue(bool_response_either.value)
        self.assertEqual(msa_source.kmers_indices.index_of("CTAGCACT"), 270)
        self.assertEqual(msa_source.kmers_indices.index_of("GACAATCA"), 1084)
        self.assertEqual(msa_source.kmers_indices.index_of("ATTTTTTT"), 1407)
        self.assertEqual(msa_source.kmers_indices.index_of("TGACAATA"), 86)
