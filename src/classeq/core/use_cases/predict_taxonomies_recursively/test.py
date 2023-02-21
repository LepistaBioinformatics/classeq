import logging
from unittest import TestCase

from classeq.core.use_cases.predict_taxonomies_recursively.calculate_joining_probability_for_group import (
    calculate_joining_probability_for_group,
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
