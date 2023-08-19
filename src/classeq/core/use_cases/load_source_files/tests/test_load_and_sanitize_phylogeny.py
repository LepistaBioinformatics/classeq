from os import getenv
from pathlib import Path
from unittest import TestCase

from classeq.core.domain.dtos.tree import ClasseqTree
from classeq.core.domain.dtos.tree_source_format import TreeSourceFormatEnum
from classeq.core.use_cases.load_source_files._load_and_sanitize_phylogeny import (
    load_and_sanitize_phylogeny,
)


class LoadAndSanitizePhylogenyTest(TestCase):
    def test_if_it_work(self) -> None:
        response = load_and_sanitize_phylogeny(
            source_file_path=Path(str(getenv("MOCK_TREE"))),
            format=TreeSourceFormatEnum.NEWICK,
            outgroups=[
                "Col_orchidophilum_BJ103_2",
                "Col_orchidophilum_CBS_119291",
                "Col_orchidophilum_IMI_309357",
                "Col_orchidophilum_BS11",
                "Col_orchidophilum_BS10",
                "Col_orchidophilum_BS06",
                "Col_orchidophilum_BS05",
                "Col_orchidophilum_CBS_63180",
                "Col_orchidophilum_CBS_63280",
                "Col_orchidophilum_COUFAL0219",
                "Col_orchidophilum_LCTJ_04",
                "Col_orchidophilum_LCTJ_02",
                "Col_orchidophilum_LCTJ_06",
                "Col_orchidophilum_LCTJ_05",
                "Col_orchidophilum_LCTJ_03",
            ],
            support_value_cutoff=70,
        )

        tree: ClasseqTree = response.value

        self.assertTrue(tree.sanitized_tree.rooted)


if __name__ == "__main__":
    from unittest import main

    main()
