import logging
from os import getenv
from pathlib import Path
from unittest import TestCase

from classeq.core.use_cases import load_and_sanitize_phylogeny

from . import train_from_single_phylogeny


class TrainFromSinglePhylogenyTest(TestCase):
    def setUp(self) -> None:
        logging.basicConfig(level=logging.DEBUG)
        self.__logger = logging.getLogger()

    def test_if_it_work(self) -> None:
        tree = load_and_sanitize_phylogeny(
            source_file_path=Path(getenv("MOCK_TREE")),  # type: ignore
            outgroups=[
                "Col_orchidophilum_CBS_119291",
                "Col_orchidophilum_IMI_309357",
                "Col_orchidophilum_CBS_63180",
                "Col_orchidophilum_CBS_63280",
            ],
        )

        train_from_single_phylogeny(tree)
