from json import dump
from pathlib import Path

import clean_base.exceptions as c_exc
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
from clean_base.either import Either, right

from classeq.core.domain.dtos.biopython_wrappers import (
    ExtendedBioPythonClade,
    ExtendedBioPythonTree,
)
from classeq.core.domain.dtos.tree import ClasseqTree
from classeq.core.domain.dtos.tree_source_format import TreeSourceFormatEnum
from classeq.settings import LOGGER, TEMP_INPUT_FILE_SUFFIX


def load_and_sanitize_phylogeny(
    source_file_path: Path,
    outgroups: list[str],
    format: TreeSourceFormatEnum,
    output_directory: Path,
    support_value_cutoff: int = 99,
) -> Either[c_exc.MappedErrors, ClasseqTree]:
    """Load phylogenetic tree into memory.

    The loading process includes the literal loading using BioPython's Phylo
    module with further reroot and sanitizing. The sanitizing process include
    the collapsing of low supported branches.

    Args:
        source_file_path (Path): The system path of the input file.
        outgroups (list[str]): A list of outgroups to reroot the phylogenetic
            tree.
        support_value_cutoff (int, optional): A cutoff for collapsing internal
            nodes of the phylogenetic tree. Defaults to 99.

    Returns:
        Either[Tree, c_exc.MappedErrors]: A phylogenetic tree in the BioPython's
            format cas the use-case was successfully executed. A default
            MappedError otherwise.
    """

    try:
        if not source_file_path.is_file():
            return c_exc.DadaTransferObjectError(
                f"Invalid path: {source_file_path}"
            )()

        # ? --------------------------------------------------------------------
        # ? Load phylogenetic tree
        #
        # Load tree and extract out-group terminals.
        #
        # ? --------------------------------------------------------------------

        LOGGER.info("Parsing and reroot tree")

        if (
            rooted_tree_either := ClasseqTree.parse_and_reroot_newick_tree(
                source_file_path,
                outgroups,
                format,
            )
        ).is_left:
            return c_exc.UseCaseError(
                "Unexpected error on parse phylogenetic tree.",
                prev=rooted_tree_either.value,
                logger=LOGGER,
            )()

        rooted_tree: Tree = rooted_tree_either.value

        # ? --------------------------------------------------------------------
        # ? Sanitize tree
        #
        # Remove nodes with low phylogenetic support.
        #
        # ? --------------------------------------------------------------------

        LOGGER.info(
            f"Collapsing low supported (< {support_value_cutoff}) branches"
        )

        sanitized_tree_either: Either = (
            ClasseqTree.collapse_low_supported_and_outgroup_nodes(
                rooted_tree=rooted_tree,
                support_value_cutoff=support_value_cutoff,
            )
        )

        if sanitized_tree_either.is_left:
            return c_exc.UseCaseError(
                "Unexpected error on sanitize phylogenetic tree.",
                prev=sanitized_tree_either.value,
                logger=LOGGER,
            )()

        sanitized_tree: Tree = sanitized_tree_either.value

        # ? --------------------------------------------------------------------
        # ? Convert default BioPython's Tree to ExtendedBioPythonTree
        # ? --------------------------------------------------------------------

        sanitized_tree = ExtendedBioPythonTree.from_bio_python_tree(
            tree=sanitized_tree,
            outgroups=outgroups,
        )

        # ? --------------------------------------------------------------------
        # ? Collapse outgroup internal nodes
        # ? --------------------------------------------------------------------

        LOGGER.info("Collapsing outgroup branches")

        tree_outgroups = [
            clade
            for clade in sanitized_tree.get_terminals()
            if clade.name in outgroups
        ]

        common_ancestor: ExtendedBioPythonClade = (
            sanitized_tree.common_ancestor(targets=tree_outgroups)
        )

        clade: ExtendedBioPythonClade
        collapsable_clades: set[ExtendedBioPythonClade] = set()

        for clade in common_ancestor.clades:
            if all([t.name in outgroups for t in clade.get_terminals()]):
                collapsable_clades.add(clade)

        sanitized_tree.collapse_all(
            lambda c: c.__hash__() in [i.__hash__() for i in collapsable_clades]
        )

        # ? --------------------------------------------------------------------
        # ? Persist sanitized tree
        # ? --------------------------------------------------------------------

        LOGGER.info("Persisting sanitized phylogenetic tree")

        cleaned_tree_file_path = output_directory.joinpath(
            "".join(
                [
                    source_file_path.stem,
                    ".",
                    TEMP_INPUT_FILE_SUFFIX,
                    ".",
                    f"cutoff{support_value_cutoff}",
                    source_file_path.suffix,
                ]
            )
        )

        LOGGER.info("Sanitized TREE as NEWICK would be persisted to:")
        LOGGER.info(f"\t{cleaned_tree_file_path.relative_to(Path.cwd())}")

        Phylo.write(
            sanitized_tree,
            cleaned_tree_file_path,
            format=format.value,
        )

        cleaned_json_file_path = cleaned_tree_file_path.with_suffix(
            "." + TreeSourceFormatEnum.PHYLO_JSON.value
        )

        LOGGER.info("Sanitized TREE as JSON would be persisted to:")
        LOGGER.info(f"\t{cleaned_json_file_path.relative_to(Path.cwd())}")

        with cleaned_json_file_path.open("w+") as f:
            dump(sanitized_tree.to_dict(), f, indent=4, default=str)

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(
            ClasseqTree(
                newick_file_path=cleaned_tree_file_path,
                phylojson_file_path=cleaned_json_file_path,
                tree_format=format,
                tree_headers=[h.name for h in sanitized_tree.get_terminals()],
                outgroups=outgroups,
                sanitized_tree=sanitized_tree,
            )
        )

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
