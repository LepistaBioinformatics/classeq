import gzip
import logging
from json import loads
from os import getenv
from pathlib import Path
from unittest import TestCase

from classeq.core.domain.dtos.reference_set import ReferenceSet


class ReferenceSetTest(TestCase):
    def setUp(self) -> None:
        logging.basicConfig(level=logging.DEBUG)
        self.__logger = logging.getLogger()

        self.__reference_set_file_path = Path(str(getenv("MOCK_REFS")))

        with gzip.open(self.__reference_set_file_path, "r") as fin:
            json_bytes = fin.read()

        json_str = json_bytes.decode("utf-8")
        reference_set_either = ReferenceSet.from_dict(content=loads(json_str))

        self.assertFalse(reference_set_either.is_left)
        self.assertTrue(reference_set_either.is_right)
        self.assertIsInstance(reference_set_either.value, ReferenceSet)

        self.__reference_set: ReferenceSet = reference_set_either.value

    def test_hierarchical_tree_generation(self) -> None:
        response_either = self.__reference_set.get_hierarchical_tree()

        self.__logger.debug(response_either)
