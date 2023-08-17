import gzip
from asyncio import run
from functools import wraps
from json import loads
from pathlib import Path
from sys import argv
from typing import Any, Callable

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
    "parse",
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
    "priors",
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
    "infer",
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
    help="The path to the priors calculated in the `priors` step.",
)
@with_resource_monitoring
def infer_identity_cmd(
    fasta_file_path: str,
    references_path: str,
    priors_path: str,
    **_: Any,
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
            )
        ).is_left:
            LOGGER.error(response.value.msg)
            click.echo("Error: Something went wrong.")
            exit(1)

    except Exception as exc:
        click.echo(f"Error: {exc}")
        exit(1)


def coro(f: Callable[..., Any]) -> Any:
    @wraps(f)
    def wrapper(*args: list[Any], **kwargs: dict[str, Any]) -> Any:
        return run(f(*args, **kwargs))

    return wrapper


""" @classeq_cmd.command(
    "serve",
    help=("Serve the phylogeny editor locally for visualize and edit TREE."),
)
@click.option(
    "-t",
    "--sanitized-tree-path",
    required=True,
    type=click.Path(
        resolve_path=True,
        readable=True,
        exists=True,
        dir_okay=True,
    ),
    help="The system path to the sanitized TREE file.",
)
@coro
async def serve_cmd(
    sanitized_tree_path: click.Path,
) -> None:
    import uvicorn
    from fastapi import FastAPI
    from fastapi.responses import FileResponse, RedirectResponse
    from fastapi.staticfiles import StaticFiles

    try:
        app = FastAPI(title="Classeq Phylogeny Editor")

        @app.get("/tree")
        async def root() -> dict[str, str]:
            return {"message": "Hello World"}

        @app.get("/")
        async def index() -> RedirectResponse:
            return RedirectResponse("/index.html")

        app.mount(
            "/",
            StaticFiles(directory="static"),
            name="static",
        )

        config = uvicorn.Config(
            app,
            port=5000,
            log_level="info",
        )

        server = uvicorn.Server(config)
        await server.serve()

    except KeyboardInterrupt:
        pass
    except Exception as exc:
        raise exc """
