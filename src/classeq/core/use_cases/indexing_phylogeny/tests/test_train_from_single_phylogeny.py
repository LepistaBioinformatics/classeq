from classeq.core.use_cases.indexing_phylogeny import (
    indexing_phylogeny,
)
from classeq.core.use_cases.indexing_phylogeny.tests.shared import (
    BaseTest,
)


class TrainFromSinglePhylogenyTest(BaseTest):
    def test_train_from_single_phylogeny(self) -> None:
        response_either = indexing_phylogeny(references=self._reference_set)

        self.assertFalse(response_either.is_left)
        self.assertTrue(response_either.is_right)


if __name__ == "__main__":
    from unittest import main

    main()
