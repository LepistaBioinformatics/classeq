from copy import deepcopy
from pathlib import Path
from typing import Any, Literal

import clean_base.exceptions as c_exc
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from clean_base.either import Either, right

from classeq.core.domain.dtos.msa_source_format import MsaSourceFormatEnum
from classeq.core.domain.dtos.priors import TreePriors
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.settings import LOGGER

from ._perform_single_sequence_phylogenetic_adherence_test import (
    perform_single_sequence_phylogenetic_adherence_test,
)


def predict_for_multiple_fasta_file(
    fasta_path: Path,
    tree_priors: TreePriors,
    reference_set: ReferenceSet,
    fasta_format: MsaSourceFormatEnum = MsaSourceFormatEnum.FASTA,
    **kwargs: Any,
) -> Either[c_exc.MappedErrors, Literal[True]]:
    print(reference_set)
    try:
        records: list[SeqRecord] = list(
            SeqIO.parse(str(fasta_path), fasta_format.value)
        )

        for record in records:
            LOGGER.debug("-" * 80)
            LOGGER.debug(f"Processing sequence: {record.id}")
            LOGGER.debug("")

            if (
                adherence_test_response_either := perform_single_sequence_phylogenetic_adherence_test(
                    target_sequence=str(record.seq),
                    reference_set=deepcopy(reference_set),
                    tree_priors=tree_priors,
                    **kwargs,
                )
            ).is_left:
                return adherence_test_response_either

            (
                adherence_final_response,
                adherence_test_path,
            ) = adherence_test_response_either.value

            LOGGER.debug("")
            LOGGER.debug(f"Adherence Terminal: {adherence_final_response}")
            LOGGER.debug("")
            LOGGER.debug(f"Adherence Test Path: {adherence_test_path}")
            LOGGER.debug("")

        return right(True)

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
