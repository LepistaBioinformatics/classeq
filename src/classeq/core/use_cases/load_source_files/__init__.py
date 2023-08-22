import gzip
from collections import defaultdict
from copy import deepcopy
from json import dump
from pathlib import Path

import clean_base.exceptions as c_exc
from attrs import asdict
from clean_base.either import Either, right

from classeq.core.domain.dtos.msa import MsaSource, MsaSourceFormatEnum
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.core.domain.dtos.tree import ClasseqTree
from classeq.core.domain.dtos.tree_source_format import TreeSourceFormatEnum
from classeq.settings import DEFAULT_KMER_SIZE, LOGGER

from ._load_and_sanitize_phylogeny import load_and_sanitize_phylogeny
from ._load_and_sanitize_sequences import load_and_sanitize_sequences


def load_source_files(
    msa_file_path: Path,
    msa_format: MsaSourceFormatEnum,
    tree_file_path: Path,
    tree_format: TreeSourceFormatEnum,
    outgroups: list[str],
    output_directory: Path | None = None,
    k_size: int = DEFAULT_KMER_SIZE,
    support_value_cutoff: int = 99,
) -> Either[c_exc.MappedErrors, ReferenceSet]:
    try:
        train_output_dir = (
            tree_file_path.parent
            if output_directory is None
            else output_directory
        )

        if not train_output_dir.is_dir():
            train_output_dir.mkdir(parents=True)

        # ? --------------------------------------------------------------------
        # ? Load MSA
        # ? --------------------------------------------------------------------

        LOGGER.info("Loading and sanitize sequences")

        if (
            msa_either := load_and_sanitize_sequences(
                source_file_path=msa_file_path,
                format=msa_format,
                output_directory=train_output_dir,
            )
        ).is_left:
            return msa_either

        msa: MsaSource = msa_either.value

        # ? --------------------------------------------------------------------
        # ? Load Tree
        # ? --------------------------------------------------------------------

        LOGGER.info("Loading and sanitize phylogeny")

        if (
            tree_either := load_and_sanitize_phylogeny(
                source_file_path=tree_file_path,
                outgroups=outgroups,
                format=tree_format,
                support_value_cutoff=support_value_cutoff,
                output_directory=train_output_dir,
            )
        ).is_left:
            return tree_either

        tree: ClasseqTree = tree_either.value

        # ? --------------------------------------------------------------------
        # ? Check content matches
        # ? --------------------------------------------------------------------

        LOGGER.info("Validating phylogeny and sequences headers")

        if not all(
            [
                *[id in msa.sequence_headers for id in tree.tree_headers],
                *[id in tree.tree_headers for id in msa.sequence_headers],
            ]
        ):
            return c_exc.UseCaseError(
                "Incompatible file contents. Tree and MSA has no "
                + "compatible identifiers.",
                exp=True,
                logger=LOGGER,
            )()

        labels_map: defaultdict[str, int] = defaultdict()

        for index, header in enumerate(msa.sequence_headers):
            labels_map[header] = index

        numeric_labels = sorted(labels_map.values())
        msa.sequence_headers = numeric_labels
        tree.tree_headers = numeric_labels

        # ? --------------------------------------------------------------------
        # ? Initialize kmers
        # ? --------------------------------------------------------------------

        LOGGER.info(f"Building kmers indices of size: {k_size}")

        if (
            init_either := msa.initialize_kmer_indices(
                headers_map=labels_map,
                k_size=k_size,
            )
        ).is_left:
            return init_either

        # ? --------------------------------------------------------------------
        # ? Persist reference set to file
        # ? --------------------------------------------------------------------

        LOGGER.info("Building output")

        references = ReferenceSet(
            tree=tree,
            msa=msa,
            labels_map=labels_map,
        )

        if (linear_tree_either := references.build_linear_tree()).is_left:
            return linear_tree_either

        tree_source = references.tree.newick_file_path

        if not train_output_dir.exists():
            train_output_dir.mkdir(parents=True)

        train_output_file_path = train_output_dir.joinpath(
            ".".join(
                [
                    tree_source.stem,
                    "reference-set",
                    f"k{k_size}",
                    "json",
                    "gz",
                ]
            )
        )

        LOGGER.info("Load output file would be persisted to:")
        LOGGER.info(f"\t{train_output_file_path}")

        with gzip.open(
            train_output_file_path, "wt", encoding="utf-8"
        ) as out_gz:
            # ? Remove sanitized tree of the persistence artifact
            tree_artifact = deepcopy(tree)
            tree_artifact.sanitized_tree = None

            dump(
                asdict(
                    ReferenceSet(
                        tree=tree_artifact,
                        msa=references.msa,
                        labels_map=references.labels_map,
                        linear_tree=references.linear_tree,
                    )
                ),
                out_gz,
                indent=4,
                default=str,
                sort_keys=True,
            )

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(references)

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
