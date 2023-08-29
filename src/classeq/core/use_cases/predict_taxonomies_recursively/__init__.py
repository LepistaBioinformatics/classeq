from collections import defaultdict
from json import dump, dumps, load
from pathlib import Path
from typing import Any, Literal
from uuid import UUID

import clean_base.exceptions as c_exc
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from clean_base.either import Either, right

from classeq.core.domain.dtos.biopython_wrappers import (
    ExtendedBioPythonClade,
    ExtendedBioPythonTree,
)
from classeq.core.domain.dtos.clade import ClasseqClade
from classeq.core.domain.dtos.msa_source_format import MsaSourceFormatEnum
from classeq.core.domain.dtos.priors import TreePriors
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.settings import LOGGER

from ._perform_adherence_test_of_child_clades import CladeAdherenceResult
from ._perform_single_sequence_phylogenetic_adherence_test import (
    perform_single_sequence_phylogenetic_adherence_test,
)


def predict_for_multiple_fasta_file(
    fasta_path: Path,
    tree_priors: TreePriors,
    reference_set: ReferenceSet,
    fasta_format: MsaSourceFormatEnum = MsaSourceFormatEnum.FASTA,
    annotated_phylojson_path: Path | None = None,
    output_file_path: Path | None = None,
    **kwargs: Any,
) -> Either[c_exc.MappedErrors, Literal[True]]:
    """Predict taxonomies for multiple sequences in a FASTA file and save the
    output as JSON and TSV.

    Args:
        fasta_path (Path): Path to the FASTA file.
        tree_priors (TreePriors): Tree priors.
        reference_set (ReferenceSet): Reference set.
        fasta_format (MsaSourceFormatEnum, optional): Format of the FASTA file.
            Defaults to MsaSourceFormatEnum.FASTA.
        annotated_phylojson_path (Path, optional): Path to the annotated
            phylojson file. Defaults to None.
        output_file_path (Path, optional): Path to the output file. Defaults to
            None.
        run_parallel (bool, optional): Whether to run the predictions in
            parallel. Defaults to True.

    returns:
        Either[c_exc.MappedErrors, Literal[True]]: Either a mapped error or
            True.

    raises:
        c_exc.UseCaseError: If an error occurs.

    """

    try:
        # ? --------------------------------------------------------------------
        # ? Resolve tree annotations
        # ? --------------------------------------------------------------------

        resolved_names: dict[UUID, ExtendedBioPythonClade] = dict()
        if annotated_phylojson_path is not None:
            with annotated_phylojson_path.open("r") as f:
                annotated_tree = ExtendedBioPythonTree.from_dict(
                    content=load(f)
                )

                resolved_names = {
                    node._id: node
                    for node in annotated_tree.get_nonterminals()
                    if node.name is not None
                }

        # ? --------------------------------------------------------------------
        # ? Predict taxonomies
        # ? --------------------------------------------------------------------

        records: list[SeqRecord] = list(
            SeqIO.parse(str(fasta_path), fasta_format.value)
        )

        response: defaultdict[str, dict[str, Any]] = defaultdict()

        if (tree_either := reference_set.get_hierarchical_tree()).is_left:
            return tree_either

        hierarchical_tree: ClasseqClade = tree_either.value

        for record in records:
            LOGGER.debug("-" * 80)
            LOGGER.debug(f"Processing sequence: {record.id}")
            LOGGER.debug("")

            if (
                adherence_test_response_either := perform_single_sequence_phylogenetic_adherence_test(
                    target_sequence=str(record.seq),
                    kmers_indices=reference_set.msa.kmers_indices,
                    kmer_size=reference_set.kmer_size,
                    tree=hierarchical_tree,
                    strand=reference_set.strand,
                    tree_priors=tree_priors,
                    **kwargs,
                )
            ).is_left:
                return adherence_test_response_either

            (
                adherence_test_path,
                status,
            ) = adherence_test_response_either.value

            resolved_path_names = __resolve_path_name(
                resolved_names=resolved_names,
                adherence_test_path=adherence_test_path,
            )

            LOGGER.debug("")
            LOGGER.debug(f"Adherence Test Path: {resolved_path_names}")
            LOGGER.debug("")

            response[record.id] = {
                "phylogeny_path": resolved_path_names,
                "status": status,
            }

            LOGGER.debug(
                "\n"
                + dumps(
                    [
                        i.to_dict(omit_children=True)
                        for i in resolved_path_names
                    ],
                    indent=4,
                    default=str,
                )
            )

        # ? --------------------------------------------------------------------
        # ? Persisting output as JSON
        # ? --------------------------------------------------------------------

        prediction_output_fields = [
            "name",
            "id",
            "support",
            "parent",
            "taxid",
            "related_rank",
        ]

        if output_file_path is None:
            output_file_path = fasta_path.parent.joinpath(
                fasta_path.stem
            ).with_suffix(".json")
        else:
            if output_file_path.suffix != ".json":
                output_file_path = output_file_path.with_suffix(".json")

        LOGGER.info(f"Saving output to `{output_file_path}`")

        if output_file_path.parent.exists() is False:
            output_file_path.parent.mkdir(parents=True)

        if output_file_path.exists():
            LOGGER.warning(f"Output file `{output_file_path}` overwritten")
            output_file_path.unlink()

        with output_file_path.open("w+") as f:
            dump(
                {
                    "input": str(fasta_path),
                    "annotations": annotated_phylojson_path,
                    "output": [
                        {
                            "target": k,
                            "status": v.get("status").value,  # type: ignore
                            "prediction": [
                                {
                                    **{
                                        y: z
                                        for y, z in i.to_dict(
                                            omit_children=True
                                        ).items()
                                        if y in prediction_output_fields
                                    },
                                    "depth": index + 1,
                                }
                                for index, i in enumerate(
                                    v.get("phylogeny_path")  # type: ignore
                                )
                            ],
                        }
                        for k, v in response.items()
                    ],
                },
                f,
                indent=4,
                default=str,
            )

        # ? --------------------------------------------------------------------
        # ? Persisting output as TSV
        # ? --------------------------------------------------------------------

        output_file_path_tsv = output_file_path.with_suffix(".tsv")

        LOGGER.info(f"Saving output to `{output_file_path_tsv}`")

        if output_file_path_tsv.exists():
            LOGGER.warning(f"Output file `{output_file_path_tsv}` overwritten")
            output_file_path_tsv.unlink()

        default_none = "NA"
        phylogeny_path: list[ClasseqClade]
        with output_file_path_tsv.open("w+") as f:
            response_lines: list[str] = []
            line_definition = "{query}\t{id}\t{status}\t{depth}\t{name}\t{taxid}\t{related_rank}\t{support}\t{branches}"
            f.write(line_definition.replace("{", "").replace("}", "") + "\n")

            for sequence_name, prediction in response.items():
                if (status := prediction.get("status")) is None:
                    continue

                if (phylogeny_path := prediction.get("phylogeny_path")) is None:  # type: ignore
                    continue

                for index, clade in enumerate(phylogeny_path):
                    response_lines.append(
                        line_definition.format(
                            query=sequence_name,
                            id=clade.id,
                            status=status.value,
                            depth=index + 1,
                            name=clade.name or default_none,
                            taxid=clade.taxid or default_none,
                            related_rank=clade.related_rank or default_none,
                            support=clade.support or default_none,
                            branches=(
                                clade.children.__len__()
                                if clade.children
                                else default_none
                            ),
                        )
                    )

            f.write("\n".join(response_lines))

        return right(True)

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()


def __resolve_path_name(
    resolved_names: dict[UUID, ExtendedBioPythonClade],
    adherence_test_path: list[CladeAdherenceResult],
) -> list[ClasseqClade]:
    renamed_clades: list[ClasseqClade] = []

    for clade_result in adherence_test_path:
        clade = clade_result.clade
        if (resolved_name := resolved_names.get(clade.id)) is not None:
            clade.name = resolved_name.name
            clade.taxid = resolved_name._taxid

            if resolved_name._related_rank is not None:
                clade.related_rank = resolved_name._related_rank.name

        renamed_clades.append(clade)
    return renamed_clades
