import gzip
from json import loads
from pathlib import Path
from sys import argv
from typing import Any

import click

from classeq.__version__ import version
from classeq.core.domain.dtos.msa_source_format import MsaSourceFormatEnum
from classeq.core.domain.dtos.priors import TreePriors
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.core.domain.dtos.tree_source_format import TreeSourceFormatEnum
from classeq.core.use_cases.load_source_files import load_source_files
from classeq.core.use_cases.predict_taxonomies_recursively import (
    predict_for_multiple_fasta_file,
)
from classeq.core.use_cases.predict_taxonomies_recursively._calculate_clade_adherence_with_bootstrap._dtos import (
    AdherenceTestStrategy,
)
from classeq.core.use_cases.train_from_single_phylogeny import (
    train_from_single_phylogeny,
)
from classeq.settings import LOGGER

from .decorators import with_resource_monitoring

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
    "std",
    help=(
        "Run standalone steps of the pipeline. Run 'classeq std --help' for "
        + "more information."
    ),
)
def std_cmd() -> None:
    pass


@classeq_cmd.group(
    "pipe",
    help=(
        "Run the complete pipeline pipeline. Run 'classeq pipe --help' for "
        + "more information."
    ),
)
def pipe_cmd() -> None:
    pass


# ? ----------------------------------------------------------------------------
# ? Initialize the CLI sub-commands
# ? ----------------------------------------------------------------------------


@std_cmd.command(
    "load",
    help=(
        "Parse a FASTA file and TREE phylogeny into a JSON file. Run "
        + "'classeq std parse --help' for more information."
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
@click.option(
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
@click.option(
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
@click.option(
    "-og",
    "--outgroups",
    type=click.STRING,
    multiple=True,
    required=False,
    help=(
        "The outgroup list. Each outgroup must be into a separate argument."
        + "Ex: -og outgroup1 -og outgroup2 -og outgroup3"
    ),
)
@click.option(
    "-of",
    "--outgroups-file",
    required=False,
    type=click.Path(
        resolve_path=True,
        readable=True,
        exists=True,
        dir_okay=True,
    ),
    help=(
        "The outgroup file path. As an alternative to the --outgroups option, "
        + "you can pass a file containing the outgroups. Each outgroup must "
        + "be into a separate line."
    ),
)
@with_resource_monitoring
def parse_source_files_cmd(
    fasta_file_path: str,
    tree_file_path: str,
    support_value_cutoff: int,
    output_directory: str | None = None,
    outgroups: tuple[str, ...] | None = None,
    outgroups_file: click.Path | None = None,
    **_: Any,
) -> None:
    """Parse a FASTA file and TREE phylogeny into a JSON file."""

    LOGGER.debug(" ".join(argv))

    if outgroups == () and outgroups_file is None:
        click.echo("Error: You must provide outgroups.")
        exit(1)

    try:
        valid_outgroups: tuple[str] = (
            outgroups
            if len(list(outgroups)) > 0  # type: ignore
            else __load_outgroups_from_file(Path(outgroups_file))  # type: ignore
        )

        if (
            response := load_source_files(
                msa_file_path=Path(fasta_file_path),
                msa_format=MsaSourceFormatEnum.FASTA,
                tree_file_path=Path(tree_file_path),
                tree_format=TreeSourceFormatEnum.NEWICK,
                outgroups=valid_outgroups,
                output_directory=(
                    Path(output_directory)
                    if output_directory is not None
                    else None
                ),
                support_value_cutoff=support_value_cutoff,
            )
        ).is_left:
            LOGGER.error(response.value.msg)
            click.echo("Error: Something went wrong.")
            exit(1)

    except Exception as exc:
        click.echo(f"Error: {exc}")
        exit(1)


def __load_outgroups_from_file(file: Path) -> list[str]:
    if not file.exists():
        raise FileNotFoundError(f"File not found: {file}")

    with open(file, "r") as f:
        outgroups = f.readlines()

    return list({outgroup.strip() for outgroup in outgroups})


@std_cmd.command(
    "train",
    help="Calculate priors of the individual phylogeny clades.",
)
@click.option(
    "-r",
    "--references",
    required=True,
    type=click.Path(
        resolve_path=True,
        readable=True,
        exists=True,
        dir_okay=True,
    ),
    help="The path to the reference file set calculated during data loading.",
)
@with_resource_monitoring
def calculate_priors_cmd(
    references: str,
    **_: Any,
) -> None:
    """Calculate priors of the individual phylogeny clades."""

    LOGGER.debug(" ".join(argv))

    try:
        references_path = Path(references)

        if not references_path.is_file():
            click.echo(
                "Error: You must provide a valid path for the reference priors."
            )

            exit(1)

        with gzip.open(references_path, "r") as fin:
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

        if (
            response := train_from_single_phylogeny(
                references=reference_set_either.value
            )
        ).is_left:
            LOGGER.error(response.value.msg)
            click.echo("Error: Something went wrong.")
            exit(1)

    except Exception as exc:
        click.echo(f"Error: {exc}")
        exit(1)


@std_cmd.command(
    "predict",
    help="Try to infer identity of multi FASTA sequences.",
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
@click.option(
    "-r",
    "--references-path",
    required=True,
    type=click.Path(
        resolve_path=True,
        readable=True,
        exists=True,
        dir_okay=True,
    ),
    help="The path to the reference file set calculated during data loading.",
)
@click.option(
    "-p",
    "--priors-path",
    required=True,
    type=click.Path(
        resolve_path=True,
        readable=True,
        exists=True,
        dir_okay=True,
    ),
    help="The path to the priors calculated in the `train` step.",
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
    "--adherence-strategy",
    required=False,
    show_default=True,
    default=AdherenceTestStrategy(None).value,
    type=click.Choice(
        [
            AdherenceTestStrategy.JOINT_PROBABILITY.value,
            AdherenceTestStrategy.KMERS_INTERSECTION.value,
        ],
        case_sensitive=False,
    ),
    help="The adherence test strategy.",
)
@click.option(
    "--calculate-bootstrap",
    is_flag=True,
    show_default=True,
    default=False,
    help="Calculate bootstrap probability of clades if True.",
)
def infer_identity_cmd(
    fasta_file_path: str,
    references_path: str,
    priors_path: str,
    adherence_strategy: str,
    annotated_phylojson_path: str | None = None,
    calculate_bootstrap: bool = False,
) -> None:
    """Try to infer identity of multi FASTA sequences."""

    LOGGER.debug(" ".join(argv))

    try:
        with gzip.open(references_path, "r") as fin:
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

        with gzip.open(priors_path, "r") as fin:
            json_bytes = fin.read()
            json_str = json_bytes.decode("utf-8")

        if (
            tree_priors_either := TreePriors.from_dict(content=loads(json_str))
        ).is_left:
            LOGGER.error(tree_priors_either.value.msg)
            click.echo("Error: Something went wrong.")
            exit(1)

        if (
            response := predict_for_multiple_fasta_file(
                fasta_path=Path(fasta_file_path),
                fasta_format=MsaSourceFormatEnum.FASTA,
                tree_priors=tree_priors_either.value,
                reference_set=reference_set_either.value,
                adherence_strategy=AdherenceTestStrategy(adherence_strategy),
                annotated_phylojson_path=(
                    Path(annotated_phylojson_path)
                    if annotated_phylojson_path is not None
                    else None
                ),
                calculate_bootstrap=calculate_bootstrap,
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
