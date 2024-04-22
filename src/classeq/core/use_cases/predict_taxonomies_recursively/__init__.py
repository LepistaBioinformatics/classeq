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
from classeq.core.domain.dtos.priors import OutgroupPriors, TreePriors
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.core.use_cases.shared.resolve_path_name import resolve_path_name
from classeq.settings import LOGGER

from ._get_outgroup_priors import get_outgroup_priors
from ._perform_single_sequence_phylogenetic_adherence_test import (
    perform_single_sequence_phylogenetic_adherence_test,
)
from ._perform_single_sequence_phylogenetic_adherence_test._dtos import (
    PredictionResult,
)


def predict_for_multiple_fasta_file(
    prediction_input_path: Path,
    tree_priors: TreePriors,
    reference_set: ReferenceSet,
    fasta_format: MsaSourceFormatEnum = MsaSourceFormatEnum.FASTA,
    outgroups_input_path: Path | None = None,
    annotated_phylojson_path: Path | None = None,
    output_file_path: Path | None = None,
    persist_as_json: bool = True,
    persist_as_csv: bool = True,
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
        # ? Calculate outgroup priors
        # ? --------------------------------------------------------------------

        outgroup_priors: OutgroupPriors | Literal[False] = False
        if outgroups_input_path is not None:
            if (
                outgroup_priors_either := get_outgroup_priors(
                    source_file_path=outgroups_input_path,
                    k_size=reference_set.kmer_size,
                    strand=reference_set.strand,
                )
            ).is_left:
                return outgroup_priors_either

            outgroup_priors = outgroup_priors_either.value

        # ? --------------------------------------------------------------------
        # ? Predict taxonomies
        # ? --------------------------------------------------------------------

        records: list[SeqRecord] = list(
            SeqIO.parse(str(prediction_input_path), fasta_format.value)
        )

        response: defaultdict[str, PredictionResult] = defaultdict()

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
                    kmer_size=reference_set.kmer_size,
                    tree=hierarchical_tree,
                    strand=reference_set.strand,
                    tree_priors=tree_priors,
                    outgroup_priors=outgroup_priors,
                    **kwargs,
                )
            ).is_left:
                return adherence_test_response_either

            prediction_result: PredictionResult = (
                adherence_test_response_either.value
            )

            resolved_path_names = resolve_path_name(
                resolved_names=resolved_names,
                adherence_test_path=prediction_result.path,
            )

            LOGGER.debug("")
            LOGGER.debug(f"Adherence Test Path: {resolved_path_names}")
            LOGGER.debug("")

            response[record.id] = prediction_result

            LOGGER.debug(
                "\n"
                + dumps(
                    [
                        i.result.to_dict(omit_children=True)
                        for i in resolved_path_names
                    ],
                    indent=4,
                    default=str,
                )
            )

        # ? --------------------------------------------------------------------
        # ? Persist output
        # ? --------------------------------------------------------------------

        if output_file_path is None:
            output_file_path = prediction_input_path.parent.joinpath(
                prediction_input_path.stem
            ).with_suffix(".json")
        else:
            if output_file_path.suffix != ".json":
                output_file_path = output_file_path.with_suffix(".json")

        if persist_as_json is True:
            __write_json(
                output_file_path=output_file_path,
                response=response,
                fasta_path=prediction_input_path,
                annotated_phylojson_path=annotated_phylojson_path,
            )

        if persist_as_csv is True:
            __write_csv(
                output_file_path=output_file_path,
                response=response,
            )

        return right(True)

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()


def __write_csv(
    output_file_path: Path,
    response: defaultdict[str, PredictionResult],
) -> None:
    output_file_path_tsv = output_file_path.with_suffix(".tsv")
    default_none = "NA"

    LOGGER.info(f"Saving output to `{output_file_path_tsv}`")

    if output_file_path_tsv.exists():
        LOGGER.warning(f"Output file `{output_file_path_tsv}` overwritten")
        output_file_path_tsv.unlink()

    with output_file_path_tsv.open("w+") as f:
        response_lines: list[str] = []
        line_definition = (
            "{query}\t{status}\t{predicted_clade_id}\t{depth}\t{name}\t{taxid}"
            + "\t{related_rank}\t{support}\t{branches}\t{query_kmers_size}"
            + "\t{in_match_kmers}\t{in_subject_kmers_size}\t{in_status}"
            + "\t{sis_match_kmers}\t{sis_subject_kmers_size}\t{sis_status}"
        )

        f.write(line_definition.replace("{", "").replace("}", "") + "\n")

        for sequence_name, prediction in response.items():
            if prediction.path.__len__() == 0:
                response_lines.append(
                    "{query}\t{status}".format(
                        query=sequence_name,
                        status=prediction.status.name,
                    )
                )

                continue

            for clade in prediction.path:
                response_lines.append(
                    line_definition.format(
                        query=sequence_name,
                        status=prediction.status.name,
                        predicted_clade_id=clade.result.clade.id,
                        depth=clade.depth,
                        name=clade.result.clade.name or default_none,
                        taxid=clade.result.clade.taxid or default_none,
                        related_rank=clade.result.clade.related_rank
                        or default_none,
                        support=clade.result.clade.support or default_none,
                        branches=(
                            clade.result.clade.children.__len__()
                            if clade.result.clade.children
                            else default_none
                        ),
                        query_kmers_size=clade.result.ingroup_adherence_test.query_kmers_size,
                        in_match_kmers=clade.result.ingroup_adherence_test.match_kmers,
                        in_subject_kmers_size=clade.result.ingroup_adherence_test.subject_kmers_size,
                        in_status=clade.result.ingroup_adherence_test.status.name,
                        sis_match_kmers=clade.result.sister_adherence_test.match_kmers,
                        sis_subject_kmers_size=clade.result.sister_adherence_test.subject_kmers_size,
                        sis_status=clade.result.sister_adherence_test.status.name,
                    )
                )

        f.write("\n".join(response_lines))


def __write_json(
    output_file_path: Path,
    response: defaultdict[str, PredictionResult],
    fasta_path: Path,
    annotated_phylojson_path: Path | None,
) -> None:
    LOGGER.info(f"Saving output to `{output_file_path}`")

    if output_file_path.parent.exists() is False:
        output_file_path.parent.mkdir(parents=True)

    if output_file_path.exists():
        LOGGER.warning(f"Output file `{output_file_path}` overwritten")
        output_file_path.unlink()

    with output_file_path.open("w+") as f:
        dump(
            {
                "input": fasta_path.name,
                "annotations": (
                    annotated_phylojson_path.name
                    if annotated_phylojson_path
                    else None
                ),
                "output": [
                    {"query": k, **v.to_dict(omit_children=True)}
                    for k, v in response.items()
                ],
            },
            f,
            indent=4,
            default=str,
        )
