from classeq.core.domain.dtos.priors import TreePriors
from classeq.core.use_cases.indexing_phylogeny._fetch_clade_specific_kmers import (
    fetch_clade_specific_kmers,
)
from classeq.core.use_cases.indexing_phylogeny.tests.shared import (
    BaseTest,
)


class EstimateCladeSpecificPriorsTest(BaseTest):
    def test_fetch_clade_specific_kmers(self) -> None:
        response_either = fetch_clade_specific_kmers(
            references=self._reference_set
        )

        self.assertFalse(response_either.is_left)
        self.assertTrue(response_either.is_right)
        self.assertIsInstance(response_either.value, TreePriors)


if __name__ == "__main__":
    from unittest import main

    main()
