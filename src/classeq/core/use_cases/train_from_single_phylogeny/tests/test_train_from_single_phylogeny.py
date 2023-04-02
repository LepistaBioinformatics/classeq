import gzip
from json import loads
from os import getenv
from pathlib import Path
from unittest import TestCase

from classeq.core.domain.dtos.priors import TreePriors
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.core.use_cases.train_from_single_phylogeny import (
    train_from_single_phylogeny,
)
from classeq.core.use_cases.train_from_single_phylogeny.estimate_clade_specific_priors import (
    estimate_clade_specific_priors,
)


class TrainFromSinglePhylogenyTest(TestCase):
    def setUp(self) -> None:
        self.__reference_set_file_path = Path(str(getenv("MOCK_REFS")))

        # ? Load reference set from mock json
        with gzip.open(self.__reference_set_file_path, "r") as fin:
            json_bytes = fin.read()

        json_str = json_bytes.decode("utf-8")
        reference_set_either = ReferenceSet.from_dict(content=loads(json_str))

        self.assertFalse(reference_set_either.is_left)
        self.assertTrue(reference_set_either.is_right)
        self.assertIsInstance(reference_set_either.value, ReferenceSet)

        self.__reference_set: ReferenceSet = reference_set_either.value

    def test_estimate_clade_specific_priors(self) -> None:
        response_either = estimate_clade_specific_priors(
            references=self.__reference_set
        )

        self.assertFalse(response_either.is_left)
        self.assertTrue(response_either.is_right)
        self.assertIsInstance(response_either.value, TreePriors)

    def test_train_from_single_phylogeny(self) -> None:
        response_either = train_from_single_phylogeny(
            references=self.__reference_set
        )

        self.assertFalse(response_either.is_left)
        self.assertTrue(response_either.is_right)


if __name__ == "__main__":
    from unittest import main

    main()
