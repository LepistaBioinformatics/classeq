from collections import defaultdict
from pathlib import Path
from uuid import uuid4

import clean_base.exceptions as c_exc
from Bio import SeqIO
from clean_base.either import Either, right

from classeq.core.domain.dtos.kmer_inverse_index import KmersInverseIndices
from classeq.core.domain.dtos.msa_source_format import MsaSourceFormatEnum
from classeq.core.domain.dtos.ordered_tuple import OrderedTuple
from classeq.core.domain.dtos.priors import (
    OutgroupLabeledPriors,
    OutgroupPriors,
)
from classeq.core.domain.dtos.strand import StrandEnum
from classeq.core.use_cases.shared.fetch_kmers_indices import (
    fetch_kmers_indices,
)
from classeq.settings import LOGGER


def get_outgroup_priors(
    source_file_path: Path,
    k_size: int,
    strand: StrandEnum,
) -> Either[c_exc.MappedErrors, OutgroupPriors]:
    try:
        # ? --------------------------------------------------------------------
        # ? Validate args
        # ? --------------------------------------------------------------------

        if not any(
            [
                source_file_path.name.endswith(suffix)
                for suffix in [".fasta", ".fna", ".fa"]
            ]
        ):
            return c_exc.UseCaseError(
                "Invalid source file. "
                + "Source file must be in .fasta, .fna, .fa format",
                logger=LOGGER,
            )()

        # ? --------------------------------------------------------------------
        # ? Generate kmer inverse indices
        # ? --------------------------------------------------------------------

        headers_map: defaultdict[str, int] = defaultdict()

        for index, record in enumerate(
            SeqIO.parse(
                handle=source_file_path,
                format=MsaSourceFormatEnum.FASTA.value,
            )
        ):
            headers_map[record.id] = -abs(index)

        if (
            kmer_indices := KmersInverseIndices.new(
                source_file_path=source_file_path,
                format=MsaSourceFormatEnum.FASTA,
                k_size=k_size,
                headers_map=headers_map,
                strand=strand,
            )
        ).is_left:
            return kmer_indices

        if (
            outgroup_kmers_either := fetch_kmers_indices(
                kmer_indices=kmer_indices.value,
                sequence_codes=[i for i in headers_map.values()],
            )
        ).is_left:
            return outgroup_kmers_either

        outgroup_priors = OutgroupPriors(
            parent=uuid4(),
            clade_priors=OutgroupLabeledPriors(
                labels=OrderedTuple([i for i in headers_map.values()]),
                kmers=outgroup_kmers_either.value,
            ),
        )

        # ? --------------------------------------------------------------------
        # ? Return a positive response
        # ? --------------------------------------------------------------------

        return right(outgroup_priors)

    except Exception as exc:
        return c_exc.UseCaseError(exc, logger=LOGGER)()
