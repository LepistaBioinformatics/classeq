from enum import Enum
from typing import Any, DefaultDict, Dict, List, Tuple
from uuid import UUID

from attrs import field, frozen


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

    def asdict(self) -> Dict[str, Any]:
        return {
            "labels": [i for i in self.labels],
            "priors": self.priors,
            "group": self.group.value,
        }


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

    def asdict(self) -> Dict[str, Any]:
        return {
            "parent": self.parent.__str__(),
            "priors": self.priors.asdict(),
        }


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

    def asdict(self) -> Dict[str, Any]:
        return {
            "parent": self.parent.__str__(),
            "priors": [i.asdict() for i in self.priors],
        }


@frozen(kw_only=True)
class TreePriors:
    """These object stores priors for the complete phylogenetic tree. Here nodes
    are refereed trough the linear tree clade IDs.
    """

    outgroup: OutgroupLabeledPriors = field()
    ingroups: List[IngroupCladePriors] = field(default=[])

    def asdict(self) -> Dict[str, Any]:
        return {
            "outgroup": self.outgroup.__str__(),
            "ingroups": [i.asdict() for i in self.ingroups],
        }
