from enum import Enum
from typing import Any, Self
from uuid import UUID

import clean_base.exceptions as c_exc
from attrs import field, frozen
from clean_base.either import Either, right

from classeq.core.domain.dtos.ordered_tuple import OrderedTuple
from classeq.settings import LOGGER


class PriorGroup(Enum):
    OUTGROUP = "outgroup"
    INGROUP = "ingroup"
    SISTER = "sister"

    def __str__(self) -> str:
        return self.name


@frozen(kw_only=True)
class CladePriors:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    labels: OrderedTuple = field()
    kmers: OrderedTuple = field()
    group: PriorGroup = field()

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return {
            "labels": self.labels,
            "kmers": self.kmers,
            "group": self.group.value,
        }

    # ? ------------------------------------------------------------------------
    # ? Public class methods
    # ? ------------------------------------------------------------------------

    @classmethod
    def from_dict(
        cls,
        content: dict[str, Any],
    ) -> Either[c_exc.MappedErrors, Self]:
        for key in ["labels", "kmers", "group"]:
            if key not in content:
                return c_exc.DadaTransferObjectError(
                    f"Invalid content detected on parse `{CladePriors}`. "
                    f"{key}` key is empty.",
                    logger=LOGGER,
                )()

        if (labels := content.get("labels")) is None:
            return c_exc.DadaTransferObjectError(
                f"Invalid content detected on parse `{CladePriors}`. "
                f"{key}` key is empty.",
                logger=LOGGER,
            )()

        if (kmers := content.get("kmers")) is None:
            return c_exc.DadaTransferObjectError(
                f"Invalid content detected on parse `{CladePriors}`. "
                f"{key}` key is empty.",
                logger=LOGGER,
            )()

        if (group := content.get("group")) is None:
            return c_exc.DadaTransferObjectError(
                f"Invalid content detected on parse `{CladePriors}`. "
                f"{key}` key is empty.",
                logger=LOGGER,
            )()

        return right(
            cls(
                labels=OrderedTuple([int(i) for i in labels]),
                kmers=OrderedTuple(kmers),
                group=PriorGroup(group.lower()),
            )
        )


@frozen(kw_only=True)
class OutgroupLabeledPriors(CladePriors):
    group: PriorGroup = field(default=PriorGroup.OUTGROUP)


@frozen(kw_only=True)
class IngroupCladePriors(CladePriors):
    group: PriorGroup = field(default=PriorGroup.INGROUP)


@frozen(kw_only=True)
class SisterCladePriors(CladePriors):
    group: PriorGroup = field(default=PriorGroup.SISTER)


@frozen(kw_only=True)
class OutgroupPriors:
    """These object stores priors of a single clade."""

    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    # The parent node denotes the most recent common ancestor (MRCA) of the
    # target clade.
    parent: UUID = field()

    # The field include the priors of all kmers shared among members of the
    # target clade.
    clade_priors: OutgroupLabeledPriors = field()

    # ? ------------------------------------------------------------------------
    # ? Validations
    # ? ------------------------------------------------------------------------

    @parent.default
    def _id_default(self) -> UUID:
        return self.default_parent_id()

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return {
            "parent": self.parent.__str__(),
            "clade_priors": self.clade_priors.to_dict(),
        }

    @staticmethod
    def default_parent_id() -> UUID:
        return UUID(int=0)

    # ? ------------------------------------------------------------------------
    # ? Public class methods
    # ? ------------------------------------------------------------------------

    @classmethod
    def from_dict(
        cls,
        content: dict[str, Any],
    ) -> Either[c_exc.MappedErrors, Self]:
        for key in [
            "parent",
            "clade_priors",
        ]:
            if key not in content:
                return c_exc.DadaTransferObjectError(
                    f"Invalid content detected on parse `{OutgroupPriors}`. "
                    f"{key}` key is empty.",
                    logger=LOGGER,
                )()

        priors_either = OutgroupLabeledPriors.from_dict(
            content=content.get("clade_priors")  # type: ignore
        )

        if priors_either.is_left:
            return priors_either

        return right(
            cls(
                parent=UUID(content.get("parent")),
                clade_priors=priors_either.value,
            )
        )


@frozen(kw_only=True)
class IngroupPriors:
    """These object stores priors of a single clade."""

    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    # The parent node denotes the most recent common ancestor (MRCA) of the
    # target clade.
    parent: UUID = field()

    # Priors of specific groups, being in order: ingroup, sister, and random
    # groups, respectively.
    clade_priors: tuple[
        IngroupCladePriors,
        SisterCladePriors,
    ] = field()

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return {
            "parent": self.parent.__str__(),
            "priors": [i.to_dict() for i in self.clade_priors],
        }

    # ? ------------------------------------------------------------------------
    # ? Public class methods
    # ? ------------------------------------------------------------------------

    @classmethod
    def from_dict(
        cls,
        content: dict[str, Any],
    ) -> Either[c_exc.MappedErrors, Self]:
        for key in [
            "parent",
            "clade_priors",
        ]:
            if key not in content:
                return c_exc.DadaTransferObjectError(
                    f"Invalid content detected on parse `{IngroupPriors}`."
                    f" {key}` key is empty.",
                    logger=LOGGER,
                )()

        priors: list[CladePriors] = []

        for prior in content.get("clade_priors"):  # type: ignore
            prior_either = CladePriors.from_dict(content=prior)

            if prior_either.is_left:
                return prior_either

            priors.append(prior_either.value)

        ingroup_prior = next(i for i in priors if i.group == PriorGroup.INGROUP)
        sister_prior = next(i for i in priors if i.group == PriorGroup.SISTER)

        return right(
            cls(
                parent=UUID(content.get("parent")),
                clade_priors=(
                    IngroupCladePriors(
                        labels=ingroup_prior.labels,
                        kmers=ingroup_prior.kmers,
                    ),
                    SisterCladePriors(
                        labels=sister_prior.labels,
                        kmers=sister_prior.kmers,
                    ),
                ),
            )
        )


@frozen(kw_only=True)
class TreePriors:
    """These object stores priors for the complete phylogenetic tree. Here nodes
    are refereed trough the linear tree clade IDs.
    """

    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    ingroups: list[IngroupPriors] = field(default=[])

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return {
            "ingroups": [i.to_dict() for i in self.ingroups],
        }

    # ? ------------------------------------------------------------------------
    # ? Public class methods
    # ? ------------------------------------------------------------------------

    @classmethod
    def from_dict(
        cls,
        content: dict[str, Any],
    ) -> Either[c_exc.MappedErrors, Self]:
        for key in ["ingroups"]:
            if key not in content:
                return c_exc.DadaTransferObjectError(
                    f"Invalid content detected on parse `{TreePriors}`. "
                    f"{key}` key is empty.",
                    logger=LOGGER,
                )()

        ingroups: list[IngroupPriors] = []
        for item in content.get("ingroups"):  # type: ignore
            parsed_item = IngroupPriors.from_dict(content=item)

            if parsed_item.is_left:
                return parsed_item

            ingroups.append(parsed_item.value)

        return right(cls(ingroups=ingroups))
