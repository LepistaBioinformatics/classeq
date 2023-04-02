from os import getenv
from pathlib import Path
from unittest import TestCase

from classeq.core.domain.dtos.msa import MsaSourceFormatEnum
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.core.domain.dtos.tree_source_format import TreeSourceFormatEnum
from classeq.core.use_cases.load_source_files import load_source_files


class LoadSourceFilesTest(TestCase):
    def test_if_it_work(self) -> None:
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
            support_value_cutoff=95,
        )

        self.assertFalse(response.is_left)
        self.assertTrue(response.is_right)
        self.assertIsInstance(response.value, ReferenceSet)


if __name__ == "__main__":
    from unittest import main

    main()
