from collections import defaultdict
from json import dump, dumps, load
from pathlib import Path
from typing import Any, Literal
from uuid import UUID

import clean_base.exceptions as c_exc
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from clean_base.either import Either, right

from classeq.core.domain.dtos.biopython_wrappers import ExtendedBioPythonTree
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
    try:
        # ? --------------------------------------------------------------------
        # ? Resolve tree annotations
        # ? --------------------------------------------------------------------

        resolved_names: dict[UUID, str] = dict()
        if annotated_phylojson_path is not None:
            with annotated_phylojson_path.open("r") as f:
                tree = ExtendedBioPythonTree.from_dict(content=load(f))

                resolved_names = {
                    node._id: node.name
                    for node in tree.get_nonterminals()
                    if node.name is not None
                }

        # ? --------------------------------------------------------------------
        # ? Predict taxonomies
        # ? --------------------------------------------------------------------

        records: list[SeqRecord] = list(
            SeqIO.parse(str(fasta_path), fasta_format.value)
        )

        response: defaultdict[str, dict[str, Any]] = defaultdict()

        for record in records:
            LOGGER.debug("-" * 80)
            LOGGER.debug(f"Processing sequence: {record.id}")
            LOGGER.debug("")

            if (
                adherence_test_response_either := perform_single_sequence_phylogenetic_adherence_test(
                    target_sequence=str(record.seq),
                    reference_set=reference_set,
                    tree_priors=tree_priors,
                    **kwargs,
                )
            ).is_left:
                return adherence_test_response_either

            (
                adherence_final_response,
                adherence_test_path,
                status,
            ) = adherence_test_response_either.value

            resolved_path_names = __resolve_path_name(
                resolved_names=resolved_names,
                adherence_test_path=adherence_test_path,
            )

            LOGGER.debug("")
            LOGGER.debug(f"Adherence Terminal: {adherence_final_response}")
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
        # ? Persisting output
        # ? --------------------------------------------------------------------

        if output_file_path is None:
            output_file_path = fasta_path.parent.joinpath(
                fasta_path.stem
            ).with_suffix(".json")
        else:
            if output_file_path.suffix != ".json":
                output_file_path = output_file_path.with_suffix(".json")

        LOGGER.info(f"Saving output to `{output_file_path}`")

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
                            "status": v.pop("status").value,
                            "prediction": [
                                {
                                    **{
                                        y: z
                                        for y, z in i.to_dict(
                                            omit_children=True
                                        ).items()
                                        if y
                                        in ["name", "id", "support", "parent"]
                                    },
                                    "depth": index + 1,
                                }
                                for index, i in enumerate(
                                    v.pop("phylogeny_path")
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

        return right(True)

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()


def __resolve_path_name(
    resolved_names: dict[UUID, str],
    adherence_test_path: list[CladeAdherenceResult],
) -> list[ClasseqClade]:
    renamed_clades: list[ClasseqClade] = []

    for clade_result in adherence_test_path:
        clade = clade_result.clade
        if (resolved_name := resolved_names.get(clade.id)) is not None:
            clade.name = resolved_name

        renamed_clades.append(clade)
    return renamed_clades
