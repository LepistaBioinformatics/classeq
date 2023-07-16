from os import getenv
from pathlib import Path
from unittest import TestCase

from classeq.core.domain.dtos.tree import TreeSource
from classeq.core.domain.dtos.tree_source_format import TreeSourceFormatEnum

from .._load_and_sanitize_phylogeny import load_and_sanitize_phylogeny


class LoadAndSanitizePhylogenyTest(TestCase):
    def test_if_it_work(self) -> None:
        response = load_and_sanitize_phylogeny(
            source_file_path=Path(str(getenv("MOCK_TREE"))),
            format=TreeSourceFormatEnum.NEWICK,
            outgroups=[
                "Col_orchidophilum_CBS_119291",
                "Col_orchidophilum_IMI_309357",
                "Col_orchidophilum_CBS_63180",
                "Col_orchidophilum_CBS_63280",
            ],
        )

        tree: TreeSource = response.value

        self.assertTrue(tree.sanitized_tree.rooted)


if __name__ == "__main__":
    from unittest import main

    main()