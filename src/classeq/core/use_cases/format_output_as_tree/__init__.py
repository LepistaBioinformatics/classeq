from json import load
from pathlib import Path
from uuid import UUID

from Bio.Phylo.NewickIO import write as write_newick
from Bio.Phylo.BaseTree import BranchColor
import clean_base.exceptions as c_exc
from clean_base.either import Either, right

from classeq.core.use_cases.predict_taxonomies_recursively._perform_single_sequence_phylogenetic_adherence_test._dtos import (
    CladeAdherenceResultStatus,
)
from classeq.core.domain.dtos.biopython_wrappers import (
    ExtendedBioPythonTree,
    ExtendedBioPythonClade,
)
from classeq.settings import (
    LOGGER,
    OUTPUT_DICT_KEY,
    PREDICTION_DICT_KEY,
    QUERY_DICT_KEY,
    STATUS_DICT_KEY,
)


def format_output_as_tree(
    phylojson_path: Path,
    results_path: Path,
    output_path: Path | None = None,
) -> Either[c_exc.MappedErrors, ...]:
    try:
        if phylojson_path is None:
            return c_exc.UseCaseError(
                "The phylojson file path is not provided", logger=LOGGER
            )()

        if results_path is None:
            return c_exc.UseCaseError(
                "The results file path is not provided", logger=LOGGER
            )()

        # ? --------------------------------------------------------------------
        # ? Load the annotated tree
        # ? --------------------------------------------------------------------

        annotated_tree: ExtendedBioPythonTree | None = None
        with phylojson_path.open("r") as f:
            annotated_tree = ExtendedBioPythonTree.from_dict(content=load(f))

        if annotated_tree is None:
            return c_exc.UseCaseError(
                "The phylojson file is empty", logger=LOGGER
            )()

        # ? --------------------------------------------------------------------
        # ? Load results
        # ? --------------------------------------------------------------------

        results = None
        with results_path.open("r") as f:
            results = load(f)

        if results is None:
            return c_exc.UseCaseError(
                "The results file is empty", logger=LOGGER
            )()

        if not isinstance(results, dict):
            return c_exc.UseCaseError(
                "The results file is not a dictionary", logger=LOGGER
            )()

        if (output_array := results.get(OUTPUT_DICT_KEY)) is None:
            return c_exc.UseCaseError(
                "The results file does not contain the output key",
                logger=LOGGER,
            )()

        if not isinstance(output_array, list):
            return c_exc.UseCaseError(
                "The output key in the results file is not a list",
                logger=LOGGER,
            )()

        # ? --------------------------------------------------------------------
        # ? Extend the tree with the results
        # ? --------------------------------------------------------------------

        for item in output_array:
            if (query := item.get(QUERY_DICT_KEY)) is None:
                return c_exc.UseCaseError(
                    f"The results file does not contain the query key: {QUERY_DICT_KEY}",
                    logger=LOGGER,
                )()

            if (status := item.get(STATUS_DICT_KEY)) is None:
                return c_exc.UseCaseError(
                    f"The results file does not contain the status key: {STATUS_DICT_KEY}",
                    logger=LOGGER,
                )()

            try:
                status = CladeAdherenceResultStatus(status)
            except ValueError:
                return c_exc.UseCaseError(
                    f"Invalid status: {status}",
                    logger=LOGGER,
                )()

            if status == CladeAdherenceResultStatus.CONCLUSIVE_OUTGROUP:
                LOGGER.warning(f"Ignored result {query} with status {status}")
                continue

            if (subject := item.get(PREDICTION_DICT_KEY)) is None:
                return c_exc.UseCaseError(
                    f"The results file does not contain the prediction key: {PREDICTION_DICT_KEY}",
                    logger=LOGGER,
                )()

            if not isinstance(subject, list):
                return c_exc.UseCaseError(
                    f"The prediction key in the results file is not a list: {subject}",
                    logger=LOGGER,
                )()

            if subject.__len__() < 1:
                LOGGER.warning(f"Ignored result {query} with empty prediction")
                continue

            terminal = subject[-1]

            if (terminal_id := terminal.get("id")) is None:
                return c_exc.UseCaseError(
                    f"The terminal does not have an id: {terminal}",
                    logger=LOGGER,
                )()

            annotated_tree.set_clade_sibling(
                id=UUID(terminal_id),
                child=ExtendedBioPythonClade(
                    name=f"QUERY__{query}",
                    branch_length=0.00001,
                    confidence=0.0,
                    color=BranchColor.from_name("black"),
                ),
            )

        # ? --------------------------------------------------------------------
        # ? Save the tree
        # ? --------------------------------------------------------------------

        if output_path is None:
            output_path = results_path.with_suffix(".newick")

        with output_path.open("w+") as f:
            write_newick([annotated_tree], f)

        return right(output_path)

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
