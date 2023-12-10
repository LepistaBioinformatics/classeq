from pathlib import Path
from typing import Any, Callable

import clean_base.exceptions as c_exc
from Bio.SeqRecord import SeqRecord
from clean_base.either import Either, right

from classeq.core.domain.dtos.clade import ClasseqClade
from classeq.core.domain.dtos.prediction_context import (
    Context,
    PredictionContext,
)
from classeq.core.domain.dtos.priors import TreePriors
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.core.use_cases.build_prediction_context.errors import (
    FailLoadingContext,
    InvalidReferenceTreeFromContext,
    UnexistingContext,
    UnexpectedErrorOnExecuteAdherenceTest,
)
from classeq.core.use_cases.predict_taxonomies_recursively._perform_single_sequence_phylogenetic_adherence_test import (
    PredictionResult,
    perform_single_sequence_phylogenetic_adherence_test,
)
from classeq.core.use_cases.shared.resolve_path_name import resolve_path_name
from classeq.settings import LOGGER


def build_prediction_context(
    context_path: Path,
) -> Either[
    c_exc.UseCaseError,
    tuple[PredictionContext, Callable[..., PredictionResult]],
]:
    """This function builds the prediction context used into a two step
    predictions scenario. The first step initialize prediction artifacts into
    memory where during the second step the prediction is performed based on the
    context delivered by the first step.
    """

    try:
        # ? --------------------------------------------------------------------
        # ? Validate arguments
        # ? --------------------------------------------------------------------

        if context_path is None:
            return c_exc.UseCaseError(
                "Invalid argument, missing 'context_path' argument.",
                logger=LOGGER,
            )()

        if not isinstance(context_path, Path):
            return c_exc.UseCaseError(
                "Invalid argument, 'context_path' argument must be a Path.",
                logger=LOGGER,
            )()

        if not context_path.exists():
            return c_exc.UseCaseError(
                f"Path {context_path} does not exist",
                logger=LOGGER,
            )()

        # ? --------------------------------------------------------------------
        # ? Load context
        # ? --------------------------------------------------------------------

        if (
            context_either := PredictionContext.load_from_yaml(
                path=context_path
            )
        ).is_left:
            return context_either
        context: PredictionContext = context_either.value

        # ? --------------------------------------------------------------------
        # ? Build callable context
        # ? --------------------------------------------------------------------

        def predict_into_context(
            model_id: str,
            seq: SeqRecord,
            **kwargs: Any,
        ) -> PredictionResult:
            named_context: Context | None
            if (
                named_context := context.get_context_by_name(
                    identifier=model_id
                )
            ) is None:
                raise UnexistingContext(
                    f"Context {model_id} does not exist.",
                )

            if named_context.prediction_priors is None:
                raise FailLoadingContext(
                    f"Context {model_id} loaded but does not have prediction priors.",
                )

            reference_set: ReferenceSet
            tree_priors: TreePriors

            (
                reference_set,
                tree_priors,
            ) = named_context.prediction_priors

            if (tree_either := reference_set.get_hierarchical_tree()).is_left:
                raise InvalidReferenceTreeFromContext(tree_either.value.msg)

            hierarchical_tree: ClasseqClade = tree_either.value

            if (
                adherence_test_response_either := perform_single_sequence_phylogenetic_adherence_test(
                    target_sequence=str(seq.seq),
                    kmer_size=reference_set.kmer_size,
                    tree=hierarchical_tree,
                    strand=reference_set.strand,
                    tree_priors=tree_priors,
                    outgroup_priors=False,
                    **kwargs,
                )
            ).is_left:
                raise UnexpectedErrorOnExecuteAdherenceTest(
                    adherence_test_response_either.value.msg
                )

            prediction_result: PredictionResult = (
                adherence_test_response_either.value
            )

            resolved_path_names = resolve_path_name(
                resolved_names=named_context.resolved_names,
                adherence_test_path=prediction_result.path,
            )

            prediction_result.path = resolved_path_names

            return prediction_result

        # ? --------------------------------------------------------------------
        # ? Return callable context
        # ? --------------------------------------------------------------------

        return right((context, predict_into_context))

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
