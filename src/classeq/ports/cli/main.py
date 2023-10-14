import gzip
import tarfile
from io import BytesIO
from json import loads
from pathlib import Path
from sys import argv
from typing import Any

import click

from classeq.__version__ import version
from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.msa_source_format import MsaSourceFormatEnum
from classeq.core.domain.dtos.priors import TreePriors
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.core.domain.dtos.strand import StrandEnum
from classeq.core.domain.dtos.tree_source_format import TreeSourceFormatEnum
from classeq.core.use_cases.indexing_phylogeny import indexing_phylogeny
from classeq.core.use_cases.load_source_files import load_source_files
from classeq.core.use_cases.predict_taxonomies_recursively import (
    predict_for_multiple_fasta_file,
)
from classeq.settings import (
    BASES,
    DEFAULT_CLASSEQ_OUTPUT_FILE_NAME,
    DEFAULT_KMER_SIZE,
    DEFAULT_MATCHES_COVERAGE,
    LOGGER,
    MINIMUM_INGROUP_QUERY_KMERS_MATCH,
    MINIMUM_INGROUP_SISTER_MATCH_KMERS_DIFFERENCE,
    REFERENCE_SET_OUTPUT_FILE_NAME,
    TRAIN_SOURCE_OUTPUT_FILE_NAME,
)

# ? ----------------------------------------------------------------------------
# ? Initialize the CLI groups
# ? ----------------------------------------------------------------------------


@click.group(
    (
        "The classeq command line interface to place sequences into real "
        + "phylogenies. Run 'classeq --help' for more information."
    )
)
@click.version_option(version)
def classeq_cmd() -> None:
    pass


@classeq_cmd.group(
    "utils",
    help="Utility commands to help you with the classeq CLI.",
)
def utils_cmd() -> None:
    pass


# ? ----------------------------------------------------------------------------
# ? Initialize the CLI sub-commands
# ? ----------------------------------------------------------------------------


__STRAND_OPTION = click.option(
    "-s",
    "--strand",
    required=False,
    type=click.Choice([s.value for s in StrandEnum]),
    default=StrandEnum(None).value,
    show_default=True,
    help="The DNA strand direction to include in analysis.",
)


__KMER_SIZE_OPTION = click.option(
    "-k",
    "--kmer-size",
    required=False,
    type=click.INT,
    default=DEFAULT_KMER_SIZE,
    show_default=True,
    help="The kmer size.",
)


__SUPPORT_VALUE_CUTOFF = click.option(
    "-s",
    "--support-value-cutoff",
    required=False,
    type=click.INT,
    default=95,
    show_default=True,
    help=(
        "The support value cutoff used to filter the phylogeny. The default "
        + "value is 95."
    ),
)


__TREE_FILE_PATH = click.option(
    "-t",
    "--tree-file-path",
    required=True,
    type=click.Path(
        resolve_path=True,
        readable=True,
        exists=True,
        dir_okay=True,
    ),
    help="The path to the TREE file.",
)


@classeq_cmd.command(
    "index",
    help=(
        "Parse a FASTA file and TREE phylogeny and calculate priors of the "
        + "individual phylogeny."
    ),
)
@click.option(
    "-f",
    "--fasta-file-path",
    required=True,
    type=click.Path(
        resolve_path=True,
        readable=True,
        exists=True,
        dir_okay=True,
    ),
    help="The path to the FASTA file.",
)
@__TREE_FILE_PATH
@click.option(
    "-o",
    "--output-directory",
    required=False,
    default=None,
    type=click.Path(
        resolve_path=True,
        readable=True,
        dir_okay=True,
    ),
    help="The path to the TREE file.",
)
@__SUPPORT_VALUE_CUTOFF
@__KMER_SIZE_OPTION
@__STRAND_OPTION
def parse_source_files_cmd(
    fasta_file_path: str,
    tree_file_path: str,
    support_value_cutoff: int,
    kmer_size: int,
    output_directory: str | None = None,
    strand: str | None = None,
    **_: Any,
) -> None:
    """Parse a FASTA file and TREE phylogeny into a JSON file."""

    LOGGER.debug(" ".join(argv))

    try:
        if (
            response := load_source_files(
                msa_file_path=Path(fasta_file_path),
                msa_format=MsaSourceFormatEnum.FASTA,
                tree_file_path=Path(tree_file_path),
                tree_format=TreeSourceFormatEnum.NEWICK,
                output_directory=(
                    Path(output_directory)
                    if output_directory is not None
                    else None
                ),
                support_value_cutoff=support_value_cutoff,
                k_size=kmer_size,
                strand=(
                    StrandEnum(None) if strand is None else StrandEnum(strand)
                ),
            )
        ).is_left:
            LOGGER.error(response.value.msg)
            click.echo("Error: Something went wrong.")
            exit(1)

        (
            indexing_file_path,
            reference_set,
        ) = response.value

        if (
            response := indexing_phylogeny(
                references=reference_set,
                indexing_file_path=indexing_file_path,
            )
        ).is_left:
            LOGGER.error(response.value.msg)
            click.echo("Error: Something went wrong.")
            exit(1)

    except Exception as exc:
        click.echo(f"Error: {exc}")
        exit(1)


@classeq_cmd.command(
    "predict",
    help="Try to infer identity of multi FASTA sequences.",
)
@click.option(
    "-q",
    "--query-fasta-path",
    required=True,
    type=click.Path(
        resolve_path=True,
        readable=True,
        exists=True,
        dir_okay=True,
    ),
    help="The path to the query FASTA file.",
)
@click.option(
    "-og",
    "--outgroups-fasta-file",
    required=False,
    type=click.Path(
        resolve_path=True,
        readable=True,
        exists=True,
        dir_okay=True,
    ),
    default=None,
    help="The outgroup file containing outgroups sequences as FASTA format.",
)
@click.option(
    "-i",
    "--classeq-indices",
    required=True,
    type=click.Path(
        resolve_path=True,
        readable=True,
        exists=True,
        dir_okay=True,
    ),
    help=f"The path to the classeq index file ({DEFAULT_CLASSEQ_OUTPUT_FILE_NAME})",
)
@click.option(
    "-o",
    "--output-file-path",
    required=False,
    type=click.Path(exists=False),
    help="The path to persist predictions output.",
)
@click.option(
    "-t",
    "--annotated-phylojson-path",
    required=False,
    type=click.Path(
        resolve_path=True,
        readable=True,
        exists=True,
        dir_okay=True,
    ),
    help="The path to the annotated phylogeny in PHYLO-JSON format.",
)
@click.option(
    "-mqm",
    "--minimum-ingroup-query-kmers-match",
    required=False,
    type=click.INT,
    default=MINIMUM_INGROUP_QUERY_KMERS_MATCH,
    show_default=True,
    help=(
        "The minimum number of matches a sequence needs to contain to be "
        + "considered belonging to a clade."
    ),
)
@click.option(
    "-isd",
    "--ingroup-sister-match-kmers-difference",
    required=False,
    type=click.INT,
    default=MINIMUM_INGROUP_SISTER_MATCH_KMERS_DIFFERENCE,
    show_default=True,
    help=(
        "The smallest difference between the ingroup and sister matches size "
        + "to consider a sequence belonging to a clade."
    ),
)
@click.option(
    "-m",
    "--matches-coverage",
    required=False,
    type=click.FLOAT,
    default=DEFAULT_MATCHES_COVERAGE,
    show_default=True,
    help=(
        "The minimum number of k-mer's matches which the query sequence needs "
        + "to contain to be considered belonging to a clade. Value must be "
        + "between 0.1 and 1."
    ),
)
def infer_identity_cmd(
    query_fasta_path: str,
    classeq_indices: str,
    outgroups_fasta_file: str | None = None,
    annotated_phylojson_path: str | None = None,
    output_file_path: str | None = None,
    **kwargs: Any,
) -> None:
    """Try to infer identity of multi FASTA sequences."""

    LOGGER.debug(" ".join(argv))

    try:
        # ? --------------------------------------------------------------------
        # ? Load the reference set
        # ? --------------------------------------------------------------------

        classeq_indices_path = Path(classeq_indices)

        if classeq_indices_path.name != DEFAULT_CLASSEQ_OUTPUT_FILE_NAME:
            raise Exception(
                f"Invalid classeq indices file name: {classeq_indices_path.name}"
            )

        with tarfile.open(classeq_indices_path) as tf:
            # ? ----------------------------------------------------------------
            # ? Load reference set
            # ? ----------------------------------------------------------------

            if (
                reference_set_content := tf.extractfile(
                    member=REFERENCE_SET_OUTPUT_FILE_NAME
                )
            ) is None:
                raise Exception(
                    f"Reference set not found: {REFERENCE_SET_OUTPUT_FILE_NAME}"
                )

            with gzip.GzipFile(
                fileobj=BytesIO(reference_set_content.read()),
                mode="rb",
            ) as fin:
                json_bytes = fin.read()
                json_str = json_bytes.decode("utf-8")

            if (
                reference_set_either := ReferenceSet.from_dict(
                    content=loads(json_str)
                )
            ).is_left:
                LOGGER.error(reference_set_either.value.msg)
                click.echo("Error: Something went wrong.")
                exit(1)

            # ? ----------------------------------------------------------------
            # ? Load indices
            # ? ----------------------------------------------------------------

            if (
                indices_content := tf.extractfile(
                    member=TRAIN_SOURCE_OUTPUT_FILE_NAME
                )
            ) is None:
                raise Exception(
                    f"Reference set not found: {TRAIN_SOURCE_OUTPUT_FILE_NAME}"
                )

            with gzip.GzipFile(
                fileobj=BytesIO(indices_content.read()),
                mode="rb",
            ) as fin:
                json_bytes = fin.read()
                json_str = json_bytes.decode("utf-8")

            if (
                tree_priors_either := TreePriors.from_dict(
                    content=loads(json_str)
                )
            ).is_left:
                LOGGER.error(tree_priors_either.value.msg)
                click.echo("Error: Something went wrong.")
                exit(1)

        # ? --------------------------------------------------------------------
        # ? Predict the taxonomies
        # ? --------------------------------------------------------------------

        if (
            response := predict_for_multiple_fasta_file(
                prediction_input_path=Path(query_fasta_path),
                outgroups_input_path=(
                    Path(outgroups_fasta_file)
                    if outgroups_fasta_file is not None
                    else outgroups_fasta_file
                ),
                fasta_format=MsaSourceFormatEnum.FASTA,
                tree_priors=tree_priors_either.value,
                reference_set=reference_set_either.value,
                annotated_phylojson_path=(
                    Path(annotated_phylojson_path)
                    if annotated_phylojson_path is not None
                    else None
                ),
                output_file_path=(
                    Path(output_file_path)
                    if output_file_path is not None
                    else None
                ),
                **kwargs,
            )
        ).is_left:
            raise Exception(response.value.msg)

    except Exception as exc:
        raise exc


@classeq_cmd.command(
    "serve",
    help=("Serve the phylogeny editor locally for visualize and edit TREE."),
)
@click.option(
    "-t",
    "--phylo-json-tree",
    required=True,
    type=click.Path(
        resolve_path=True,
        readable=True,
        exists=True,
        dir_okay=True,
    ),
    help="The system path to the sanitized TREE file as PHYLO-JSON format.",
)
def serve_cmd(
    phylo_json_tree: str,
) -> None:
    from classeq.ports.app.main import main

    main(phylo_json_tree=Path(phylo_json_tree))


@utils_cmd.command(
    "get-kmers",
    help=(
        "Get kmers of a character sequence. Only DNA/RNA residuals should be "
        + f"accepted ({', '.join(BASES)})."
    ),
)
@__STRAND_OPTION
@__KMER_SIZE_OPTION
def get_sequence_kmers_cmd(
    strand: str,
    kmer_size: int,
) -> None:
    sequence = click.get_text_stream("stdin").read()

    for line in sequence.splitlines():
        if line.startswith(">"):
            click.get_text_stream("stdout").write(f"{line}\n")
            continue

        for kmer in KmersInverseIndices.generate_kmers(
            dna_sequence=KmersInverseIndices.sanitize_sequence(
                sequence=line,
                as_seq=True,
            ),
            k_size=kmer_size or DEFAULT_KMER_SIZE,
            strand=StrandEnum(strand),
        ):
            click.get_text_stream("stdout").write(f"{kmer}\n")


@utils_cmd.command(
    "sanitize-tree",
    help="Reroot tree at midpoint and remove low supported branches.",
)
@__TREE_FILE_PATH
@__SUPPORT_VALUE_CUTOFF
@click.option(
    "-o",
    "--output-file",
    required=True,
    type=click.Path(
        resolve_path=True,
        readable=True,
        dir_okay=True,
    ),
    help="The file name to persist tree.",
)
def sanitize_tree_cmd(
    tree_file_path: str,
    support_value_cutoff: int,
    output_file: str,
) -> None:
    from Bio.Phylo.BaseTree import Tree
    from Bio import Phylo

    from classeq.core.domain.dtos.biopython_wrappers import (
        ExtendedBioPythonTree,
    )
    from classeq.core.domain.dtos.tree import ClasseqTree

    if (
        rooted_tree_either := ClasseqTree.parse_and_reroot_newick_tree(
            Path(tree_file_path),
            TreeSourceFormatEnum.NEWICK,
        )
    ).is_left:
        raise Exception(rooted_tree_either.value)

    rooted_tree: Tree = rooted_tree_either.value

    rooted_tree_extended = ExtendedBioPythonTree.from_bio_python_tree(
        tree=rooted_tree,
    )

    sanitized_tree_either = ClasseqTree.collapse_low_supported_nodes(
        rooted_tree=rooted_tree_extended,
        support_value_cutoff=support_value_cutoff,
    )

    if sanitized_tree_either.is_left:
        raise Exception(sanitized_tree_either.value)

    sanitized_tree: Tree = sanitized_tree_either.value

    Phylo.write(
        sanitized_tree,
        Path(output_file),
        format=TreeSourceFormatEnum.NEWICK.value,
    )
