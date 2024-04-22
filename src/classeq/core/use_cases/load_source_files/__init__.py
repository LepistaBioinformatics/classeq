import gzip
import tarfile
from collections import defaultdict
from json import dump
from pathlib import Path
from tempfile import TemporaryDirectory

import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.domain.dtos.msa import MsaSource, MsaSourceFormatEnum
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.core.domain.dtos.strand import StrandEnum
from classeq.core.domain.dtos.tree import ClasseqTree
from classeq.core.domain.dtos.tree_source_format import TreeSourceFormatEnum
from classeq.settings import (
    DEFAULT_CLASSEQ_OUTPUT_FILE_NAME,
    LOGGER,
    REFERENCE_SET_OUTPUT_FILE_NAME,
)

from ._load_and_sanitize_phylogeny import load_and_sanitize_phylogeny
from ._load_and_sanitize_sequences import load_and_sanitize_sequences


def load_source_files(
    msa_file_path: Path,
    msa_format: MsaSourceFormatEnum,
    tree_file_path: Path,
    tree_format: TreeSourceFormatEnum,
    k_size: int,
    strand: StrandEnum,
    output_directory: Path | None = None,
    support_value_cutoff: int = 99,
    rescale_to_100: bool = False,
) -> Either[c_exc.MappedErrors, tuple[Path, ReferenceSet]]:
    try:
        train_output_dir = (
            tree_file_path.parent
            if output_directory is None
            else output_directory
        )

        if not train_output_dir.is_dir():
            train_output_dir.mkdir(parents=True)

        train_output_dir_output_file = train_output_dir.joinpath(
            DEFAULT_CLASSEQ_OUTPUT_FILE_NAME
        )

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

        msa: MsaSource
        sanitized_msa_path: Path

        (
            sanitized_msa_path,
            msa,
        ) = msa_either.value

        # ? --------------------------------------------------------------------
        # ? Load Tree
        # ? --------------------------------------------------------------------

        LOGGER.info("Loading and sanitize phylogeny")

        if (
            tree_either := load_and_sanitize_phylogeny(
                source_file_path=tree_file_path,
                format=tree_format,
                support_value_cutoff=support_value_cutoff,
                output_directory=train_output_dir,
                rescale_to_100=rescale_to_100,
            )
        ).is_left:
            return tree_either

        classeq_tree: ClasseqTree = tree_either.value

        # ? --------------------------------------------------------------------
        # ? Check content matches
        # ? --------------------------------------------------------------------

        LOGGER.info("Validating phylogeny and sequences headers")

        if not all(
            [
                *[
                    id in msa.sequence_headers
                    for id in classeq_tree.tree_headers
                ],
                *[
                    id in classeq_tree.tree_headers
                    for id in msa.sequence_headers
                ],
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
        classeq_tree.tree_headers = numeric_labels

        # ? --------------------------------------------------------------------
        # ? Initialize kmers
        # ? --------------------------------------------------------------------

        LOGGER.info(f"Building kmers indices of size: {k_size}")

        if (
            init_either := msa.initialize_kmer_indices(
                headers_map=labels_map,
                k_size=k_size,
                source_file_path=sanitized_msa_path,
                strand=strand,
            )
        ).is_left:
            return init_either

        # ? --------------------------------------------------------------------
        # ? Persist reference set to file
        # ? --------------------------------------------------------------------

        LOGGER.info("Building output")

        classeq_tree.update_tree_hash()

        references = ReferenceSet(
            kmer_size=k_size,
            strand=strand,
            tree=classeq_tree,
            msa=msa,
            labels_map=labels_map,
        )

        if (linear_tree_either := references.build_linear_tree()).is_left:
            return linear_tree_either

        if not train_output_dir.exists():
            train_output_dir.mkdir(parents=True)

        temp_dir = TemporaryDirectory()

        train_output_file_path = Path(temp_dir.name).joinpath(
            REFERENCE_SET_OUTPUT_FILE_NAME
        )

        with gzip.open(train_output_file_path, "wt", encoding="utf-8") as out:
            dump(
                references.to_dict(),
                out,
                indent=4,
                default=str,
                sort_keys=True,
            )

        with tarfile.open(train_output_dir_output_file, "w:") as tar:
            tar.add(
                train_output_file_path,
                arcname=train_output_file_path.name,
            )

        temp_dir.cleanup()

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right((train_output_dir_output_file, references))

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
