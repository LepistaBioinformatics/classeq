import gzip
from json import loads
from os import getenv
from pathlib import Path
from unittest import TestCase

from classeq.core.domain.dtos.reference_set import ReferenceSet


class BaseTest(TestCase):
    def setUp(self) -> None:
        reference_set_file_path = Path(str(getenv("MOCK_REFS")))

        # ? Load reference set from mock json
        with gzip.open(reference_set_file_path, "r") as fin:
            json_bytes = fin.read()

        json_str = json_bytes.decode("utf-8")
        reference_set_either = ReferenceSet.from_dict(content=loads(json_str))

        self.assertFalse(reference_set_either.is_left)
        self.assertTrue(reference_set_either.is_right)
        self.assertIsInstance(reference_set_either.value, ReferenceSet)

        self._reference_set: ReferenceSet = reference_set_either.value
