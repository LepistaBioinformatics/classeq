from hashlib import md5
from json import load
from pathlib import Path
from typing import Any, Self

import clean_base.exceptions as c_exc
from attr import define, field
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
from clean_base.either import Either, right
from classeq.core.domain.dtos.biopython_wrappers import ExtendedBioPythonTree

from classeq.core.domain.dtos.msa_source_format import MsaSourceFormatEnum
from classeq.core.domain.dtos.tree_source_format import TreeSourceFormatEnum
from classeq.settings import LOGGER


@define(kw_only=True)
class ClasseqTree:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    tree_headers: list[int] = field()
    tree_format: MsaSourceFormatEnum = field()
    sanitized_tree: ExtendedBioPythonTree = field()
    tree_hash: str | None = field(default=None)

    # ? ------------------------------------------------------------------------
    # ? Public class methods
    # ? ------------------------------------------------------------------------

    @classmethod
    def from_dict(
        cls,
        content: dict[str, Any],
    ) -> Either[c_exc.MappedErrors, Self]:
        for key in [
            "tree_headers",
            "tree_format",
            "sanitized_tree",
            "tree_hash",
        ]:
            if key not in content:
                return c_exc.DadaTransferObjectError(
                    f"Invalid content detected on parse `{ClasseqTree}`. "
                    f"{key}` key is empty.",
                    logger=LOGGER,
                )()

        if (sanitized_tree := content.get("sanitized_tree")) is None:
            return c_exc.DadaTransferObjectError(
                f"Invalid content detected on parse `{ClasseqTree}`. "
                f"`sanitized_tree` key is empty.",
                logger=LOGGER,
            )()

        tree = cls(
            tree_headers=content.get("tree_headers"),  # type: ignore
            tree_format=eval(content.get("tree_format")),  # type: ignore
            sanitized_tree=ExtendedBioPythonTree.from_dict(
                content=sanitized_tree,
            ),
        )

        if (tree_hash := content.get("tree_hash")) is None:
            return c_exc.DadaTransferObjectError(
                f"Invalid content detected on parse `{ClasseqTree}`. "
                f"`tree_hash` key is empty.",
                logger=LOGGER,
            )()

        tree.update_tree_hash()

        if tree.tree_hash != tree_hash:
            return c_exc.DadaTransferObjectError(
                f"Invalid content detected on parse `{ClasseqTree}`. "
                f"`tree_hash` key is empty.",
                logger=LOGGER,
            )()

        return right(tree)

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return {
            "tree_headers": self.tree_headers,
            "tree_format": self.tree_format,
            "sanitized_tree": self.sanitized_tree.to_dict(),
            "tree_hash": self.tree_hash,
        }

    def update_tree_hash(self) -> None:
        self.tree_hash = md5(
            str(
                (
                    self.tree_headers,
                    self.tree_format,
                    "".join(
                        [
                            f"{clade._id}{clade.branch_length}{clade.confidence}"
                            for clade in self.sanitized_tree.root.get_nonterminals()
                        ]
                    ),
                    "".join(
                        [
                            f"{clade._id}{clade.branch_length}{clade.name}"
                            for clade in self.sanitized_tree.root.get_terminals()
                        ]
                    ),
                )
            ).encode()
        ).hexdigest()

    # ? ------------------------------------------------------------------------
    # ? Public static methods
    # ? ------------------------------------------------------------------------

    @staticmethod
    def parse_and_reroot_phylojson_tree(
        phylojson_file_path: Path,
    ) -> Either[c_exc.MappedErrors, Tree]:
        try:
            with phylojson_file_path.open("r") as phylojson_file:
                content = load(phylojson_file)
                raw_tree: ExtendedBioPythonTree = (
                    ExtendedBioPythonTree.from_dict(content)
                )

            raw_tree.root_at_midpoint()

            return right(raw_tree)

        except Exception as exc:
            return c_exc.UseCaseError(exc, logger=LOGGER)()

    @staticmethod
    def parse_and_reroot_newick_tree(
        newick_file_path: Path,
        format: TreeSourceFormatEnum,
    ) -> Either[c_exc.MappedErrors, Tree]:
        try:
            raw_tree: Tree = Phylo.read(newick_file_path, format.value)
            raw_tree.root_at_midpoint()

            return right(raw_tree)

        except Exception as exc:
            return c_exc.UseCaseError(exc, logger=LOGGER)()

    @staticmethod
    def collapse_low_supported_nodes(
        rooted_tree: Tree,
        support_value_cutoff: int = 99,
    ) -> Either[c_exc.MappedErrors, Tree]:
        try:
            if rooted_tree.rooted is False:
                return c_exc.DadaTransferObjectError(
                    "Un-rooted trees not allowed."
                )()

            rooted_tree.collapse_all(
                lambda c: c.confidence is not None
                and c.confidence < support_value_cutoff
            )

            return right(rooted_tree)

        except Exception as exc:
            return c_exc.UseCaseError(exc, logger=LOGGER)()
