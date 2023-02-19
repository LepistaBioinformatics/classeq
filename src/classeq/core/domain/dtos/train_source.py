from attr import define, field

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.dtos.msa import MsaSource
from classeq.core.domain.dtos.tree import TreeSource
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER


@define
class TrainSource:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    tree: TreeSource = field()
    fasta: MsaSource = field()

    # ? ------------------------------------------------------------------------
    # ? Validations
    # ? ------------------------------------------------------------------------

    # TODO: do implement fasta and phylogeny validations

    # ? ------------------------------------------------------------------------
    # ? Public Instance Methods
    # ? ------------------------------------------------------------------------

    def check_attrs_compatibility(self) -> Either[bool, c_exc.MappedErrors]:
        try:
            return right(False)

        except Exception as exc:
            return left(c_exc.UseCaseError(exc, logger=LOGGER))
