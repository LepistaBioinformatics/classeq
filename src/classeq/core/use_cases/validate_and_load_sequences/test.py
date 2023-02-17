import logging
from os import getenv
from pathlib import Path
from unittest import TestCase

from classeq.core.domain.dtos.msa import MsaSourceFormatEnum

from . import validate_and_load_sequences


class LoadSequencesAndCountKmerTest(TestCase):
    def setUp(self) -> None:
        logging.basicConfig(level=logging.DEBUG)
        self.__logger = logging.getLogger()

    def test_if_it_work(self) -> None:
        response = validate_and_load_sequences(
            source_file_path=Path(getenv("MOCK_MSA")),
            format=MsaSourceFormatEnum.FASTA,
        )

        self.__logger.debug(response)
