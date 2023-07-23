from os import getenv
from pathlib import Path
from unittest import TestCase

from classeq.core.domain.dtos.msa import MsaSource, MsaSourceFormatEnum
from classeq.core.use_cases.load_source_files._load_and_sanitize_sequences import (
    load_and_sanitize_sequences,
)


class LoadAndSanitizeSequencesTest(TestCase):
    def test_if_it_work(self) -> None:
        response = load_and_sanitize_sequences(
            source_file_path=Path(str(getenv("MOCK_MSA"))),
            format=MsaSourceFormatEnum.FASTA,
        )

        self.assertFalse(response.is_left)
        self.assertTrue(response.is_right)
        self.assertIsInstance(response.value, MsaSource)


if __name__ == "__main__":
    from unittest import main

    main()
