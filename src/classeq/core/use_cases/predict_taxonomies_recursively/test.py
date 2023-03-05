import gzip
import logging
from json import loads
from os import getenv
from pathlib import Path
from unittest import TestCase

from classeq.core.domain.dtos.priors import TreePriors
from classeq.core.use_cases.predict_taxonomies_recursively.calculate_joining_probability_for_group import (
    calculate_joining_probability_for_group,
)
from classeq.core.use_cases.predict_taxonomies_recursively.do_clade_adherence_test_for_single_sequence import (
    do_clade_adherence_test_for_single_sequence,
)
from classeq.core.use_cases.predict_taxonomies_recursively.perform_phylogenetic_adherence_test import (
    perform_phylogenetic_adherence_test,
)


class PredictTaxonomiesRecursivelyTest(TestCase):
    def setUp(self) -> None:
        logging.basicConfig(level=logging.DEBUG)
        self.__logger = logging.getLogger()

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

        self.__Col_cuscutae_CSL_473 = (
            "CCTTCATTGAGACCAAGTACGCTGTGAGTATCACCCCCACTTTACCCCTCCATAATGATATCACGTCTGCTACAATAACACCAGCTTCATCGGTAACCACGGGAAAAGAGTCAGAGCTAGTACTCTCGACTCTTTGGACCCAAGGTTTCGATTGGGCTCGTTGTTGTAATGATACGACGTGACACAATCATGCAGAAACAGCCCAAACAAAATTTGCTGACAGACAATCATCACAGGCCTACATGCTCAAGTAC",
        )

        self.__Col_orchidophilum_IMI_309357 = (
            "CCTTCATTGAGACCAAGTACGCTGTGAGTATCACCCCACTTTACCCCTCCATGATGATATCACATCTGTCACGACAATACCAGCCTCATCGGCCACTGGGAAAGAAATGAGCTAGCACTCTCGATCCTGTGACCCAGGATACTGAAGCGGCTCGTCCCAATGGCATGATGTGACTAGGTCACGAAGAAATAGTTGGGACAACATTTGCTGACAGACCACTACCACAGGCCTACATGCTCAAGTAC",
        )

        self.__Col_acutatum_CBS_110735 = (
            "CCTTCATTGAGACCAAGTACGCTGTGAGTATCACCCCACTTTACCCCTCCATCATGATATCACGTCTGCCACGATAACACCAGCTTCGTCGGTACCCACGGGAAAAGAGTCAGAGCTAGCGCTCTCGACTCTTTTGCCCCGAGGTTTCGATTGGGCTCGTTGTAATGATGCGACGTGATACAACCATGCAGAAACAGCCGAGACAAAATTTGCTGACAGACAATCATCACAGGCCTACATGCTCAAGTAC",
        )

        self.__priors_file_path = Path(str(getenv("MOCK_PRIORS")))
        self.__reference_set_file_path = Path(str(getenv("MOCK_REFS")))

    def test_calculate_joining_probability_for_group_work(self) -> None:
        response = calculate_joining_probability_for_group(
            self.__prediction_target,
            self.__ingroup_reference,
            self.__sister_reference,
        )

        self.__logger.debug(response.value)

    def test_perform_phylogenetic_adherence_test_work(self) -> None:
        response = perform_phylogenetic_adherence_test(
            self.__prediction_target,
            self.__ingroup_reference,
            self.__sister_reference,
            self.__outgroup_reference,
        )

        self.__logger.debug(response.value)

    def test_do_clade_adherence_test_for_single_sequence(self) -> None:
        with gzip.open(self.__priors_file_path, "r") as fin:
            json_bytes = fin.read()

        json_str = json_bytes.decode("utf-8")
        tree_priors_either = TreePriors.from_dict(content=loads(json_str))

        self.assertFalse(tree_priors_either.is_left)
        self.assertTrue(tree_priors_either.is_right)
        self.assertIsInstance(tree_priors_either.value, TreePriors)

        do_clade_adherence_test_for_single_sequence(
            target_sequence=self.__Col_cuscutae_CSL_473,
            clade_priors="",
            kmer_indices="",
        )
