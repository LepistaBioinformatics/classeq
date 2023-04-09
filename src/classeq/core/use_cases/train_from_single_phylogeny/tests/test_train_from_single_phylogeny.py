from classeq.core.use_cases.train_from_single_phylogeny import (
    train_from_single_phylogeny,
)
from classeq.core.use_cases.train_from_single_phylogeny.tests.shared import (
    BaseTest,
)


class TrainFromSinglePhylogenyTest(BaseTest):
    def test_train_from_single_phylogeny(self) -> None:
        response_either = train_from_single_phylogeny(
            references=self._reference_set
        )

        self.assertFalse(response_either.is_left)
        self.assertTrue(response_either.is_right)


if __name__ == "__main__":
    from unittest import main

    main()
