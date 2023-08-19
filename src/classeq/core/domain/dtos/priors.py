from collections import defaultdict
from enum import Enum
from typing import Any, Self
from uuid import UUID

import clean_base.exceptions as c_exc
from attrs import field, frozen
from clean_base.either import Either, left, right

from classeq.settings import LOGGER


class PriorGroup(Enum):
    OUTGROUP = "outgroup"
    INGROUP = "ingroup"
    SISTER = "sister"
    NOISE = "noise"

    def __str__(self) -> str:
        return self.name


@frozen(kw_only=True)
class LabeledPriors:
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    labels: tuple[int, ...] = field()
    priors: defaultdict[str, float] = field()
    group: PriorGroup = field()

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return {
            "labels": [i for i in self.labels],
            "priors": self.priors,
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
        for key in [
            "labels",
            "priors",
            "group",
        ]:
            if key not in content:
                return left(
                    c_exc.DadaTransferObjectError(
                        f"Invalid content detected on parse `{LabeledPriors}`. "
                        f"{key}` key is empty.",
                        logger=LOGGER,
                    )
                )

        return right(
            cls(
                labels=tuple([int(i) for i in content.get("labels")]),  # type: ignore
                priors=defaultdict(float, content.get("priors")),  # type: ignore
                group=PriorGroup(content.get("group").lower()),  # type: ignore
            )
        )


@frozen(kw_only=True)
class OutgroupLabeledPriors(LabeledPriors):
    group: PriorGroup = field(default=PriorGroup.OUTGROUP)


@frozen(kw_only=True)
class IngroupLabeledPriors(LabeledPriors):
    group: PriorGroup = field(default=PriorGroup.INGROUP)


@frozen(kw_only=True)
class SisterGroupLabeledPriors(LabeledPriors):
    group: PriorGroup = field(default=PriorGroup.SISTER)


@frozen(kw_only=True)
class NoiseGroupLabeledPriors(LabeledPriors):
    group: PriorGroup = field(default=PriorGroup.NOISE)


@frozen(kw_only=True)
class OutgroupCladePriors:
    """These object stores priors of a single clade."""

    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    # The parent node denotes the most recent common ancestor (MRCA) of the
    # target clade.
    parent: UUID = field()

    # The field include the priors of all kmers shared among members of the
    # target clade.
    priors: OutgroupLabeledPriors = field()

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return {
            "parent": self.parent.__str__(),
            "priors": self.priors.to_dict(),
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
            "priors",
        ]:
            if key not in content:
                return left(
                    c_exc.DadaTransferObjectError(
                        f"Invalid content detected on parse `{OutgroupCladePriors}`. "
                        f"{key}` key is empty.",
                        logger=LOGGER,
                    )
                )

        priors_either = OutgroupLabeledPriors.from_dict(
            content=content.get("priors")  # type: ignore
        )

        if priors_either.is_left:
            return priors_either

        return right(
            cls(
                parent=UUID(content.get("parent")),
                priors=priors_either.value,
            )
        )


@frozen(kw_only=True)
class IngroupCladePriors:
    """These object stores priors of a single clade."""

    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    # The parent node denotes the most recent common ancestor (MRCA) of the
    # target clade.
    parent: UUID = field()

    # Priors of specific groups, being in order: ingroup, sister, and random
    # groups, respectively.
    priors: tuple[
        IngroupLabeledPriors,
        SisterGroupLabeledPriors,
    ] = field()

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return {
            "parent": self.parent.__str__(),
            "priors": [i.to_dict() for i in self.priors],
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
            "priors",
        ]:
            if key not in content:
                return left(
                    c_exc.DadaTransferObjectError(
                        f"Invalid content detected on parse `{IngroupCladePriors}`."
                        f" {key}` key is empty.",
                        logger=LOGGER,
                    )
                )

        priors: list[LabeledPriors] = []

        for prior in content.get("priors"):  # type: ignore
            prior_either = LabeledPriors.from_dict(content=prior)

            if prior_either.is_left:
                return prior_either

            priors.append(prior_either.value)

        ingroup_prior = next(i for i in priors if i.group == PriorGroup.INGROUP)
        sister_prior = next(i for i in priors if i.group == PriorGroup.SISTER)

        return right(
            cls(
                parent=UUID(content.get("parent")),
                priors=(
                    IngroupLabeledPriors(
                        labels=ingroup_prior.labels,
                        priors=ingroup_prior.priors,
                    ),
                    SisterGroupLabeledPriors(
                        labels=sister_prior.labels,
                        priors=sister_prior.priors,
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

    outgroup: OutgroupCladePriors = field()
    ingroups: list[IngroupCladePriors] = field(default=[])

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return {
            "outgroup": self.outgroup.to_dict(),
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
        for key in [
            "outgroup",
            "ingroups",
        ]:
            if key not in content:
                return left(
                    c_exc.DadaTransferObjectError(
                        f"Invalid content detected on parse `{TreePriors}`. "
                        f"{key}` key is empty.",
                        logger=LOGGER,
                    )
                )

        outgroup_either = OutgroupCladePriors.from_dict(
            content=content.get("outgroup")  # type: ignore
        )

        if outgroup_either.is_left:
            return outgroup_either

        ingroups: list[IngroupCladePriors] = []
        for item in content.get("ingroups"):  # type: ignore
            parsed_item = IngroupCladePriors.from_dict(content=item)

            if parsed_item.is_left:
                return parsed_item

            ingroups.append(parsed_item.value)

        return right(
            cls(
                outgroup=outgroup_either.value,
                ingroups=ingroups,
            )
        )
