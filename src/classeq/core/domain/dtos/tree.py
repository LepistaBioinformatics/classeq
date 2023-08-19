from json import load
from pathlib import Path
from typing import Any, Self

import clean_base.exceptions as c_exc
from attr import define, field
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
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

    newick_file_path: Path = field()
    phylojson_file_path: Path = field()
    tree_headers: list[int] = field()
    outgroups: list[str] = field()
    tree_format: MsaSourceFormatEnum = field()
    sanitized_tree: ExtendedBioPythonTree | None = field(default=None)

    # ? ------------------------------------------------------------------------
    # ? Public class methods
    # ? ------------------------------------------------------------------------

    @classmethod
    def from_dict(
        cls,
        content: dict[str, Any],
    ) -> Either[c_exc.MappedErrors, Self]:
        for key in [
            "newick_file_path",
            "phylojson_file_path",
            "tree_headers",
            "outgroups",
            "tree_format",
            "sanitized_tree",
        ]:
            if key not in content:
                return c_exc.DadaTransferObjectError(
                    f"Invalid content detected on parse `{ClasseqTree}`. "
                    f"{key}` key is empty.",
                    logger=LOGGER,
                )()

        return right(
            cls(
                newick_file_path=Path(content.get("newick_file_path")),  # type: ignore
                phylojson_file_path=Path(content.get("phylojson_file_path")),  # type: ignore
                tree_headers=content.get("tree_headers"),  # type: ignore
                outgroups=content.get("outgroups"),  # type: ignore
                tree_format=eval(content.get("tree_format")),  # type: ignore
                sanitized_tree=content.get("sanitized_tree"),  # type: ignore
            )
        )

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def load_tree(
        self,
        outgroups: list[str],
    ) -> Either[c_exc.MappedErrors, bool]:
        try:
            if not self.newick_file_path.is_file():
                return c_exc.DadaTransferObjectError(
                    f"Invalid path: {self.newick_file_path}"
                )()

            # ? ----------------------------------------------------------------
            # ? Load phylogenetic tree
            # ? ----------------------------------------------------------------

            LOGGER.info("Loading phylogenetic tree")

            if (
                rooted_tree_either := self.parse_and_reroot_phylojson_tree(
                    self.phylojson_file_path,
                    outgroups,
                )
            ).is_left:
                return c_exc.UseCaseError(
                    "Unexpected error on parse phylogenetic tree.",
                    prev=rooted_tree_either.value,
                    logger=LOGGER,
                )()

            self.sanitized_tree = rooted_tree_either.value

            return right(True)

        except Exception as exc:
            return c_exc.UseCaseError(exc, logger=LOGGER)()

    # ? ------------------------------------------------------------------------
    # ? Public static methods
    # ? ------------------------------------------------------------------------

    @staticmethod
    def parse_and_reroot_phylojson_tree(
        phylojson_file_path: Path,
        outgroups: list[str],
    ) -> Either[c_exc.MappedErrors, Tree]:
        try:
            with phylojson_file_path.open("r") as phylojson_file:
                content = load(phylojson_file)
                raw_tree: ExtendedBioPythonTree = (
                    ExtendedBioPythonTree.from_dict(content)
                )

            terminals: list[Clade] = raw_tree.get_terminals()

            tree_outgroups = [
                clade for clade in terminals if clade.name in outgroups
            ]

            outgroup_paths: list[Clade] = list()
            for clade in tree_outgroups:
                ancestors: list[Clade] = raw_tree.get_path(clade)

                outgroup_paths.extend(
                    [i for i in ancestors if i.confidence is not None]
                )

            raw_tree.collapse_all(lambda c: c in outgroup_paths)

            if not all([out.name in outgroups for out in tree_outgroups]):
                return c_exc.UseCaseError(
                    f"Not all specified outgroups exists in phylogeny: {outgroups}",
                    exp=True,
                    logger=LOGGER,
                )()

            if raw_tree.root_with_outgroup(tree_outgroups) is None:
                LOGGER.warning("Outgroup is the current tree root")
                raw_tree.rooted = True

            return right(raw_tree)

        except Exception as exc:
            return c_exc.UseCaseError(exc, logger=LOGGER)()

    @staticmethod
    def parse_and_reroot_newick_tree(
        newick_file_path: Path,
        outgroups: list[str],
        format: TreeSourceFormatEnum,
    ) -> Either[c_exc.MappedErrors, Tree]:
        try:
            raw_tree: Tree = Phylo.read(newick_file_path, format.value)

            terminals: list[Clade] = raw_tree.get_terminals()

            tree_outgroups = [
                clade for clade in terminals if clade.name in outgroups
            ]

            outgroup_paths: list[Clade] = list()
            for clade in tree_outgroups:
                ancestors: list[Clade] = raw_tree.get_path(clade)

                outgroup_paths.extend(
                    [i for i in ancestors if i.confidence is not None]
                )

            raw_tree.collapse_all(lambda c: c in outgroup_paths)

            if not all([out.name in outgroups for out in tree_outgroups]):
                return c_exc.UseCaseError(
                    f"Not all specified outgroups exists in phylogeny: {outgroups}",
                    exp=True,
                    logger=LOGGER,
                )()

            if raw_tree.root_with_outgroup(tree_outgroups) is None:
                LOGGER.warning("Outgroup is the current tree root")
                raw_tree.rooted = True

            return right(raw_tree)

        except Exception as exc:
            return c_exc.UseCaseError(exc, logger=LOGGER)()

    @staticmethod
    def collapse_low_supported_and_outgroup_nodes(
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
