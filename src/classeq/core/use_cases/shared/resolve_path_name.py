from uuid import UUID

from classeq.core.domain.dtos.biopython_wrappers import (
    ExtendedBioPythonClade,
)
from classeq.core.use_cases.predict_taxonomies_recursively._perform_single_sequence_phylogenetic_adherence_test._dtos import (
    PredictionStep,
)


def resolve_path_name(
    resolved_names: dict[UUID, ExtendedBioPythonClade],
    adherence_test_path: list[PredictionStep],
) -> list[PredictionStep]:
    renamed_clades = []

    for clade_result in adherence_test_path:
        clade = clade_result.result.clade
        if (resolved_name := resolved_names.get(clade.id)) is not None:
            clade.name = resolved_name.name
            clade.taxid = resolved_name._taxid

            if resolved_name._related_rank is not None:
                clade.related_rank = resolved_name._related_rank.name

        renamed_clades.append(clade_result)
    return renamed_clades
