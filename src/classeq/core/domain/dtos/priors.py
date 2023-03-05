from collections import defaultdict
from enum import Enum
from typing import Any, DefaultDict, Dict, List, Self, Tuple
from uuid import UUID

from attrs import field, frozen

import classeq.core.domain.utils.exceptions as c_exc
from classeq.core.domain.utils.either import Either, left, right
from classeq.settings import LOGGER


class PriorGroup(Enum):
    OUTGROUP = "outgroup"
    INGROUP = "ingroup"
    SISTER = "sister"
    NOISE = "noise"


@frozen(kw_only=True)
class LabeledPriors:
    labels: Tuple[int, ...] = field()
    priors: DefaultDict[str, float] = field()
    group: PriorGroup = field()

    def to_dict(self) -> Dict[str, Any]:
        return {
            "labels": [i for i in self.labels],
            "priors": self.priors,
            "group": self.group.value,
        }

    @classmethod
    def from_dict(
        cls,
        content: Dict[str, Any],
    ) -> Either[Self, c_exc.MappedErrors]:
        for key in [
            "labels",
            "priors",
            "group",
        ]:
            if key not in content:
                return left(
                    c_exc.InvalidArgumentError(
                        f"Invalid content detected on parse `{LabeledPriors}`. "
                        f"{key}` key is empty.",
                        logger=LOGGER,
                    )
                )

        return right(
            cls(
                labels=tuple([int(i) for i in content.get("labels")]),  # type: ignore
                priors=defaultdict(float, content.get("priors")),  # type: ignore
                group=PriorGroup(content.get("group")),
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
class CladePriors:
    """These object stores priors of a single clade."""

    # The parent node denotes the most recent common ancestor (MRCA) of the
    # target clade.
    parent: UUID = field()

    # The field include the priors of all kmers shared among members of the
    # target clade.
    priors: LabeledPriors = field()

    def to_dict(self) -> Dict[str, Any]:
        return {
            "parent": self.parent.__str__(),
            "priors": self.priors.to_dict(),
        }

    @classmethod
    def from_dict(
        cls,
        content: Dict[str, Any],
    ) -> Either[Self, c_exc.MappedErrors]:
        for key in [
            "parent",
            "priors",
        ]:
            if key not in content:
                return left(
                    c_exc.InvalidArgumentError(
                        f"Invalid content detected on parse `{CladePriors}`. "
                        f"{key}` key is empty.",
                        logger=LOGGER,
                    )
                )

        priors_either = LabeledPriors.from_dict(
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

    # The parent node denotes the most recent common ancestor (MRCA) of the
    # target clade.
    parent: UUID = field()

    # Priors of specific groups, being in order: ingroup, sister, and random
    # groups, respectively.
    priors: Tuple[
        IngroupLabeledPriors,
        SisterGroupLabeledPriors,
        NoiseGroupLabeledPriors,
    ] = field()

    def to_dict(self) -> Dict[str, Any]:
        return {
            "parent": self.parent.__str__(),
            "priors": [i.to_dict() for i in self.priors],
        }

    @classmethod
    def from_dict(
        cls,
        content: Dict[str, Any],
    ) -> Either[Self, c_exc.MappedErrors]:
        for key in [
            "parent",
            "priors",
        ]:
            if key not in content:
                return left(
                    c_exc.InvalidArgumentError(
                        f"Invalid content detected on parse `{IngroupCladePriors}`."
                        f" {key}` key is empty.",
                        logger=LOGGER,
                    )
                )

        priors: List[LabeledPriors] = []

        for prior in content.get("priors"):  # type: ignore
            prior_either = LabeledPriors.from_dict(content=prior)

            if prior_either.is_left:
                return prior_either

            priors.append(prior_either.value)

        ingroup_prior = next(i for i in priors if i.group == PriorGroup.INGROUP)
        sister_prior = next(i for i in priors if i.group == PriorGroup.SISTER)
        noise_prior = next(i for i in priors if i.group == PriorGroup.NOISE)

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
                    NoiseGroupLabeledPriors(
                        labels=noise_prior.labels,
                        priors=noise_prior.priors,
                    ),
                ),
            )
        )


@frozen(kw_only=True)
class TreePriors:
    """These object stores priors for the complete phylogenetic tree. Here nodes
    are refereed trough the linear tree clade IDs.
    """

    outgroup: CladePriors = field()
    ingroups: List[IngroupCladePriors] = field(default=[])

    def to_dict(self) -> Dict[str, Any]:
        return {
            "outgroup": self.outgroup.to_dict(),
            "ingroups": [i.to_dict() for i in self.ingroups],
        }

    @classmethod
    def from_dict(
        cls,
        content: Dict[str, Any],
    ) -> Either[Self, c_exc.MappedErrors]:
        for key in [
            "outgroup",
            "ingroups",
        ]:
            if key not in content:
                return left(
                    c_exc.InvalidArgumentError(
                        f"Invalid content detected on parse `{TreePriors}`. "
                        f"{key}` key is empty.",
                        logger=LOGGER,
                    )
                )

        outgroup_either = CladePriors.from_dict(
            content=content.get("outgroup")  # type: ignore
        )

        if outgroup_either.is_left:
            return outgroup_either

        ingroups: List[IngroupCladePriors] = []
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
