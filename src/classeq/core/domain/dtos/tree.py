from pathlib import Path
from typing import List

from attr import define, field
from Bio.Phylo.BaseTree import Tree


@define(kw_only=True)
class TreeSource:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    source_file_path: Path = field()
    tree_headers: List[str] = field()
    raw_sequences: Tree | None = field(init=False, default=None)

    # ? ------------------------------------------------------------------------
    # ? Public class methods
    # ? ------------------------------------------------------------------------

    # TODO: do implement `new` classmethod
