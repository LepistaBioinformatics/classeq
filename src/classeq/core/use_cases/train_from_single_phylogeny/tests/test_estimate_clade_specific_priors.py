from classeq.core.domain.dtos.priors import TreePriors
from classeq.core.use_cases.train_from_single_phylogeny._estimate_clade_specific_priors import (
    estimate_clade_specific_priors,
)
from classeq.core.use_cases.train_from_single_phylogeny.tests.shared import (
    BaseTest,
)


class EstimateCladeSpecificPriorsTest(BaseTest):
    def test_estimate_clade_specific_priors(self) -> None:
        response_either = estimate_clade_specific_priors(
            references=self._reference_set
        )

        self.assertFalse(response_either.is_left)
        self.assertTrue(response_either.is_right)
        self.assertIsInstance(response_either.value, TreePriors)


if __name__ == "__main__":
    from unittest import main

    main()
