import gzip
from json import loads
from os import getenv
from pathlib import Path
from unittest import TestCase

from classeq.core.domain.dtos.priors import TreePriors
from classeq.core.domain.dtos.reference_set import ReferenceSet


class BaseTest(TestCase):
    def setUp(self) -> None:
        self._Col_cuscutae_CSL_473 = "CCTTCATTGAGACCAAGTACGCTGTGAGTATCACCCCCACTTTACCCCTCCATAATGATATCACGTCTGCTACAATAACACCAGCTTCATCGGTAACCACGGGAAAAGAGTCAGAGCTAGTACTCTCGACTCTTTGGACCCAAGGTTTCGATTGGGCTCGTTGTTGTAATGATACGACGTGACACAATCATGCAGAAACAGCCCAAACAAAATTTGCTGACAGACAATCATCACAGGCCTACATGCTCAAGTAC"
        self._Col_orchidophilum_IMI_309357 = "CCTTCATTGAGACCAAGTACGCTGTGAGTATCACCCCACTTTACCCCTCCATGATGATATCACATCTGTCACGACAATACCAGCCTCATCGGCCACTGGGAAAGAAATGAGCTAGCACTCTCGATCCTGTGACCCAGGATACTGAAGCGGCTCGTCCCAATGGCATGATGTGACTAGGTCACGAAGAAATAGTTGGGACAACATTTGCTGACAGACCACTACCACAGGCCTACATGCTCAAGTAC"
        self._Col_acutatum_CBS_110735 = "CCTTCATTGAGACCAAGTACGCTGTGAGTATCACCCCACTTTACCCCTCCATCATGATATCACGTCTGCCACGATAACACCAGCTTCGTCGGTACCCACGGGAAAAGAGTCAGAGCTAGCGCTCTCGACTCTTTTGCCCCGAGGTTTCGATTGGGCTCGTTGTAATGATGCGACGTGATACAACCATGCAGAAACAGCCGAGACAAAATTTGCTGACAGACAATCATCACAGGCCTACATGCTCAAGTAC"
        self._Col_paranaense_CBS_134729 = "CCTTCATTGAGACCAAGTACGCTGTGAGTATCACCCCACTTTACCCCTCCATCATGATATCACGTCTGCCACGATAACACCAGCTTCGTCGGTACCCACGGCAAAAGAGTCAGGACTAGCACTCTCGACTTTTTTGCCCCAGGGTTTCGATTGGGCTTGTTGTAATGACACGACGTGACACAATCATGCAGAAACAGCCGAGACAGAACTTGCTGACAGACAATCATCACAGGCCTACATGCTCAAGTAC"
        self._Col_nymphaeae_CBS_100064 = "ccttcattgagaccaagtacgctgtgagtatcaccccactttacccctccatcatgatatcacgtctgccacgataacaccagcttcgtcgatatccacgggaaaagagtcggagctagcactctcgactcttttgtcccaaggtttcgattgggcttgttgtaacgacacgacgtgacacaatcatgcagaaacagccgagacaaaacttgctgacagacaatcatcacaggcctacatgctcaagtac"
        self._Col_godetiae_CBS_16050 = "ccttcattgagaccaagtacgctgtgagtatcaccccactttacccctccatgatgatatcacgtctgtcacgataacaccaccctaatcggtaaccatgggaaagagccagagctgctagcactctcgactcttttcccccaaggtttagatttggctcgttgcaatggcaagacgtgacgagatcatgtagaaacatccaagacaaaatttgctgacagacaatcatcacaggcctacatgctcaagtac"
        self._Col_salicis_CBS_19156 = "CCTTCATTGAGACCAAGTACGCTGTTAGTATCACCCCACTTTACCCCCCCCCCAATGATGATATCACGTCTGCCACGTTAACACCACCCTAATCGGTAACCACGGGAAAGAGCCAGAGCTGCTAGCACTCTCGACTCTTTTGCCCCAAGGTTTCGATTTGGCTCGTTGCAATTGGCACGACGTGATGGGATCATGTAGAAACACCCAAGACAATATTTGCTGACAGACAATCATCACAGGCCTACATGCTCAAGTAC"
        self._Col_simmondsii_CBS_111531 = "CCTTCATTGAGACCAAGTACGCTGTGAGTATCACCCCACTTTACCCCTCCATCATGCTATCACGTCTACCACGATAACACCAGCTTCGTCGTTATCCGCGGGGAAAAGAGTGAGAGCTAGCAATCTCGACTCTTTGACCCCAAGGTTTCGATTGGGCTTGTTGTAATGACACGACGTGACACAATCGTGCAGAAACAGCCGAGACCAAACTTGCTGACAGACAATTCATCACAGGCCTACATGCTCAAGTAC"

        self._priors_file_path = Path(str(getenv("MOCK_PRIORS")))
        self._reference_set_file_path = Path(str(getenv("MOCK_REFS")))

        # ? Load tree priors from mock json
        with gzip.open(self._priors_file_path, "r") as fin:
            json_bytes = fin.read()

        json_str = json_bytes.decode("utf-8")
        tree_priors_either = TreePriors.from_dict(content=loads(json_str))

        self.assertFalse(tree_priors_either.is_left)
        self.assertTrue(tree_priors_either.is_right)
        self.assertIsInstance(tree_priors_either.value, TreePriors)

        self._tree_priors: TreePriors = tree_priors_either.value

        # ? Load reference set from mock json
        with gzip.open(self._reference_set_file_path, "r") as fin:
            json_bytes = fin.read()

        json_str = json_bytes.decode("utf-8")
        reference_set_either = ReferenceSet.from_dict(content=loads(json_str))

        self.assertFalse(reference_set_either.is_left)
        self.assertTrue(reference_set_either.is_right)
        self.assertIsInstance(reference_set_either.value, ReferenceSet)

        self._reference_set: ReferenceSet = reference_set_either.value