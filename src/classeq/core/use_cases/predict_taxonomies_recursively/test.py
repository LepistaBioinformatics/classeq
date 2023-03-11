import gzip
from json import loads
from os import getenv
from pathlib import Path
from unittest import TestCase

from classeq.core.domain.dtos.priors import TreePriors
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.core.use_cases.predict_taxonomies_recursively.do_clade_adherence_test_for_single_sequence import (
    do_clade_adherence_test_for_single_sequence,
)


class PredictTaxonomiesRecursivelyTest(TestCase):
    def setUp(self) -> None:
        self.__ingroup_reference = [
            ["AAAA", "AAAC", "AAAT", "AAAG"],
            ["AAAA", "AAAT", "AAAG", "AACA"],
            ["AAAA", "AAAT", "AAAG", "AACA"],
            ["AAAA", "TAAT", "AAAG", "AACA"],
            ["AAAA", "TAAT", "AAAG", "AACA"],
            ["AAAA", "AAAT", "AAAG", "AACA"],
            ["AAAA", "AAAT", "AAAG", "AACA"],
            ["AAAA", "AAAT", "AAAG", "AACA"],
            ["AAAA", "AAAT", "AAAG", "AACA"],
        ]

        self.__sister_reference = [
            ["CAAA", "AAAC", "GAAT", "AAAG"],
            ["CAAA", "AAAT", "GAAG", "AACA"],
            ["CAAA", "AAAT", "GAAG", "AACA"],
            ["AAAA", "TAAT", "AAAG", "AACA"],
            ["AAAA", "TAAT", "AAAG", "AACA"],
            ["AAAA", "TAAT", "AAAG", "AACA"],
            ["AAAA", "TAAT", "AAAG", "AACA"],
            ["AAAA", "TAAT", "AAAG", "AACA"],
            ["AAAA", "TAAT", "AAAG", "AACA"],
            ["AAAA", "TAAT", "AAAG", "AACA"],
        ]

        self.__outgroup_reference = [
            ["CAAA", "AAAC", "GAAT", "AAAG"],
            ["CAAA", "AAAT", "GAAG", "AACA"],
            ["CAAA", "AAAC", "GAAT", "AAAG"],
            ["CAAA", "AAAC", "GAAT", "AAAG"],
            ["CAAA", "AAAT", "GAAG", "AACA"],
            ["AAAA", "TAAT", "CAAG", "TACA"],
            ["AAAA", "TAAT", "CAAG", "TACA"],
            ["AAAA", "TAAT", "CAAG", "TACA"],
            ["AAAA", "TAAT", "CAAG", "TACA"],
            ["AAAA", "TAAT", "CAAG", "TACA"],
        ]

        self.__prediction_target = [
            "AAAA",
            "AAAT",
            "AACA",
            "AACA",
        ]

        self.__Col_cuscutae_CSL_473 = "CCTTCATTGAGACCAAGTACGCTGTGAGTATCACCCCCACTTTACCCCTCCATAATGATATCACGTCTGCTACAATAACACCAGCTTCATCGGTAACCACGGGAAAAGAGTCAGAGCTAGTACTCTCGACTCTTTGGACCCAAGGTTTCGATTGGGCTCGTTGTTGTAATGATACGACGTGACACAATCATGCAGAAACAGCCCAAACAAAATTTGCTGACAGACAATCATCACAGGCCTACATGCTCAAGTAC"
        self.__Col_orchidophilum_IMI_309357 = "CCTTCATTGAGACCAAGTACGCTGTGAGTATCACCCCACTTTACCCCTCCATGATGATATCACATCTGTCACGACAATACCAGCCTCATCGGCCACTGGGAAAGAAATGAGCTAGCACTCTCGATCCTGTGACCCAGGATACTGAAGCGGCTCGTCCCAATGGCATGATGTGACTAGGTCACGAAGAAATAGTTGGGACAACATTTGCTGACAGACCACTACCACAGGCCTACATGCTCAAGTAC"
        self.__Col_acutatum_CBS_110735 = "CCTTCATTGAGACCAAGTACGCTGTGAGTATCACCCCACTTTACCCCTCCATCATGATATCACGTCTGCCACGATAACACCAGCTTCGTCGGTACCCACGGGAAAAGAGTCAGAGCTAGCGCTCTCGACTCTTTTGCCCCGAGGTTTCGATTGGGCTCGTTGTAATGATGCGACGTGATACAACCATGCAGAAACAGCCGAGACAAAATTTGCTGACAGACAATCATCACAGGCCTACATGCTCAAGTAC"

        self.__priors_file_path = Path(str(getenv("MOCK_PRIORS")))
        self.__reference_set_file_path = Path(str(getenv("MOCK_REFS")))

        # ? Load tree priors from mock json
        with gzip.open(self.__priors_file_path, "r") as fin:
            json_bytes = fin.read()

        json_str = json_bytes.decode("utf-8")
        tree_priors_either = TreePriors.from_dict(content=loads(json_str))

        self.assertFalse(tree_priors_either.is_left)
        self.assertTrue(tree_priors_either.is_right)
        self.assertIsInstance(tree_priors_either.value, TreePriors)

        self.__tree_priors: TreePriors = tree_priors_either.value

        # ? Load reference set from mock json
        with gzip.open(self.__reference_set_file_path, "r") as fin:
            json_bytes = fin.read()

        json_str = json_bytes.decode("utf-8")
        reference_set_either = ReferenceSet.from_dict(content=loads(json_str))

        self.assertFalse(reference_set_either.is_left)
        self.assertTrue(reference_set_either.is_right)
        self.assertIsInstance(reference_set_either.value, ReferenceSet)

        self.__reference_set: ReferenceSet = reference_set_either.value

    def test_do_clade_adherence_test_for_single_sequence(self) -> None:
        outgroup_clade_adherence_either = (
            do_clade_adherence_test_for_single_sequence(
                target_sequence=self.__Col_cuscutae_CSL_473,
                clade_priors=self.__tree_priors.outgroup,
                kmer_indices=self.__reference_set.msa.kmers_indices,
            )
        )

        self.assertFalse(outgroup_clade_adherence_either.is_left)
        self.assertTrue(outgroup_clade_adherence_either.is_right)

        for prior in self.__tree_priors.ingroups:
            response = do_clade_adherence_test_for_single_sequence(
                target_sequence=self.__Col_cuscutae_CSL_473,
                clade_priors=prior,
                kmer_indices=self.__reference_set.msa.kmers_indices,
            )

            self.assertFalse(response.is_left)
            self.assertTrue(response.is_right)
