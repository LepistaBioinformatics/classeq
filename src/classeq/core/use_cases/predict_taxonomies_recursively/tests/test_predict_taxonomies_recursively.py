from classeq.core.use_cases.predict_taxonomies_recursively._do_clade_adherence_test_for_single_sequence import (
    do_clade_adherence_test_for_single_sequence,
)
from classeq.core.use_cases.predict_taxonomies_recursively.tests.shared import (
    BaseTest,
)


class PredictTaxonomiesRecursivelyTest(BaseTest):
    def test_do_clade_adherence_test_for_single_sequence(self) -> None:
        outgroup_clade_adherence_either = (
            do_clade_adherence_test_for_single_sequence(
                target_sequence=self._Col_cuscutae_CSL_473,
                clade_priors=self._tree_priors.outgroup,
                kmer_indices=self._reference_set.msa.kmers_indices,
                total_length=len(self._reference_set.labels_map),
            )
        )

        self.assertFalse(outgroup_clade_adherence_either.is_left)
        self.assertTrue(outgroup_clade_adherence_either.is_right)

        for prior in self._tree_priors.ingroups:
            response = do_clade_adherence_test_for_single_sequence(
                target_sequence=self._Col_cuscutae_CSL_473,
                clade_priors=prior,
                kmer_indices=self._reference_set.msa.kmers_indices,
                total_length=len(self._reference_set.labels_map),
            )

            self.assertFalse(response.is_left)
            self.assertTrue(response.is_right)


if __name__ == "__main__":
    from unittest import main

    main()
