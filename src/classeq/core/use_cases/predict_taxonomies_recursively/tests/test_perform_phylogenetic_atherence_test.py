from classeq.core.use_cases.predict_taxonomies_recursively._perform_phylogenetic_adherence_test import (
    perform_phylogenetic_adherence_test,
)
from classeq.core.use_cases.predict_taxonomies_recursively.tests.shared import (
    BaseTest,
)


class PerformPhylogeneticAdherenceTest(BaseTest):
    def test_if_works(self) -> None:
        outgroup = [
            (
                "Col_orchidophilum_IMI_309357",
                self._Col_orchidophilum_IMI_309357,
            )
        ]

        largest_clade = [
            # ("Col_acutatum_CBS_110735", self._Col_acutatum_CBS_110735),
            ("Col_simmondsii_CBS_111531", self._Col_simmondsii_CBS_111531),
            # ("Col_paranaense_CBS_134729", self._Col_paranaense_CBS_134729),
            # ("Col_nymphaeae_CBS_100064", self._Col_nymphaeae_CBS_100064),
            # ("Col_cuscutae_CSL_473", self._Col_cuscutae_CSL_473),
            # ("Col_godetiae_CBS_16050", self._Col_godetiae_CBS_16050),
        ]

        minor_clade = [
            ("Col_salicis_CBS_19156", self._Col_salicis_CBS_19156),
        ]

        for name, seq in [
            *outgroup,
            *largest_clade,
            *minor_clade,
        ]:
            print(f"\n{name}")
            response = perform_phylogenetic_adherence_test(
                target_sequence=seq,
                reference_set=self._reference_set,
                tree_priors=self._tree_priors,
            )

            self.assertFalse(response.is_left)
            self.assertTrue(response.is_right)


if __name__ == "__main__":
    from unittest import main

    main()
