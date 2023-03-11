from copy import deepcopy
import gzip
from collections import defaultdict
from json import dump
from pathlib import Path
from typing import DefaultDict

from attrs import asdict

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.msa import MsaSource, MsaSourceFormatEnum
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.core.domain.dtos.tree import TreeSource
from classeq.core.domain.dtos.tree_source_format import TreeSourceFormatEnum
from classeq.core.domain.utils.either import Either, left, right
from classeq.core.use_cases.load_source_files.load_and_sanitize_phylogeny import (
    load_and_sanitize_phylogeny,
)
from classeq.core.use_cases.load_source_files.load_and_sanitize_sequences import (
    load_and_sanitize_sequences,
)
from classeq.settings import LOGGER


def load_source_files(
    msa_file_path: Path,
    msa_format: MsaSourceFormatEnum,
    tree_file_path: Path,
    tree_format: TreeSourceFormatEnum,
    outgroups: list[str],
    support_value_cutoff: int = 99,
) -> Either[bool, c_exc.MappedErrors]:
    try:
        # ? --------------------------------------------------------------------
        # ? Load MSA
        # ? --------------------------------------------------------------------

        msa_either = load_and_sanitize_sequences(
            source_file_path=msa_file_path,
            format=msa_format,
        )

        if msa_either.is_left:
            return msa_either

        msa: MsaSource = msa_either.value

        # ? --------------------------------------------------------------------
        # ? Load Tree
        # ? --------------------------------------------------------------------

        tree_either = load_and_sanitize_phylogeny(
            source_file_path=tree_file_path,
            outgroups=outgroups,
            format=tree_format,
            support_value_cutoff=support_value_cutoff,
        )

        if tree_either.is_left:
            return tree_either

        tree: TreeSource = tree_either.value

        # ? --------------------------------------------------------------------
        # ? Check content matches
        # ? --------------------------------------------------------------------

        if not all(
            [
                *[id in msa.sequence_headers for id in tree.tree_headers],
                *[id in tree.tree_headers for id in msa.sequence_headers],
            ]
        ):
            return left(
                c_exc.InvalidArgumentError(
                    "Incompatible file contents. Tree and MSA has no "
                    + "compatible identifiers.",
                    exp=True,
                    logger=LOGGER,
                )
            )

        # ? --------------------------------------------------------------------
        # ? Check content matches
        # ? --------------------------------------------------------------------

        labels_map: DefaultDict[str, int] = defaultdict()

        for index, header in enumerate(msa.sequence_headers):
            labels_map[header] = index

        numeric_labels = sorted(labels_map.values())
        msa.sequence_headers = numeric_labels
        tree.tree_headers = numeric_labels

        # ? --------------------------------------------------------------------
        # ? Initialize kmers
        # ? --------------------------------------------------------------------

        init_either = msa.initialize_kmer_indices(headers_map=labels_map)

        if init_either.is_left:
            return init_either

        # ? --------------------------------------------------------------------
        # ? Persist reference set to file
        # ? --------------------------------------------------------------------

        references = ReferenceSet(
            tree=tree,
            msa=msa,
            labels_map=labels_map,
        )

        tree_source = references.tree.source_file_path
        train_output_file_path = tree_source.parent.joinpath(
            ".".join(
                [
                    tree_source.stem,
                    "reference-set",
                    "json",
                    "gz",
                ]
            )
        )

        LOGGER.info(
            f"Train output file would be persisted to: {train_output_file_path}"
        )

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
                    )
                ),
                out_gz,
                indent=4,
                default=str,
                sort_keys=True,
            )

        linear_tree_output_file_path = tree_source.parent.joinpath(
            ".".join(
                [
                    tree_source.stem,
                    "linear-tree",
                    "json",
                    "gz",
                ]
            )
        )

        LOGGER.info(
            f"Linear tree would be persisted to: {linear_tree_output_file_path}"
        )

        with gzip.open(
            linear_tree_output_file_path, "wt", encoding="utf-8"
        ) as out_gz:
            linear_tree = references.get_linear_tree()

            if linear_tree.is_left:
                return linear_tree

            dump(
                [asdict(i) for i in linear_tree.value],
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
        return left(c_exc.UseCaseError(exc, logger=LOGGER))
