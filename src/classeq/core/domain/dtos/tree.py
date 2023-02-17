from enum import Enum
from pathlib import Path
from typing import List, Optional

from attr import define, field
from Bio.Phylo.BaseTree import Tree


class TreeSourceFormatEnum(Enum):
    NEWICK = "newick"


@define(kw_only=True)
class TreeSource:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    source_file_path: Path = field()
    tree_headers: List[str] = field()
    raw_sequences: Optional[Tree] = field()

    # ? ------------------------------------------------------------------------
    # ? Public class methods
    # ? ------------------------------------------------------------------------

    # TODO: do implement `new`` classmethod
