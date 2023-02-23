import logging
from os import getenv
from pathlib import Path
from unittest import TestCase

from classeq.core.domain.dtos.msa_source_format import MsaSourceFormatEnum
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.core.domain.dtos.tree_source_format import TreeSourceFormatEnum
from classeq.core.use_cases.load_source_files import load_source_files
from classeq.core.use_cases.train_from_single_phylogeny.estimate_clade_specific_conditional_probabilities import (
    estimate_clade_specific_conditional_probabilities,
)
from classeq.core.use_cases.train_from_single_phylogeny.estimate_global_kmer_specific_priors import (
    estimate_global_kmer_specific_priors,
)


class TrainFromSinglePhylogenyTest(TestCase):
    def setUp(self) -> None:
        logging.basicConfig(level=logging.DEBUG)
        self.__logger = logging.getLogger()

        msa_file_path = Path(str(getenv("MOCK_MSA")))
        tree_file_path = Path(str(getenv("MOCK_TREE")))

        response = load_source_files(
            msa_file_path=msa_file_path,
            msa_format=MsaSourceFormatEnum.FASTA,
            tree_file_path=tree_file_path,
            tree_format=TreeSourceFormatEnum.NEWICK,
            outgroups=[
                "Col_orchidophilum_CBS_119291",
                "Col_orchidophilum_IMI_309357",
                "Col_orchidophilum_CBS_63180",
                "Col_orchidophilum_CBS_63280",
            ],
            support_value_cutoff=70,
        )

        self.__ref: ReferenceSet = response.value

    def test_estimate_global_kmer_specific_priors(self) -> None:
        response_either = estimate_global_kmer_specific_priors(
            msa=self.__ref.msa,
        )

        self.assertFalse(response_either.is_left)
        self.assertTrue(response_either.is_right)

    def test_estimate_clade_specific_conditional_probabilities(self) -> None:
        estimate_clade_specific_conditional_probabilities(references=self.__ref)
