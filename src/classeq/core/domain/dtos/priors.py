from attrs import frozen, field
from uuid import UUID
from typing import DefaultDict, List


@frozen(kw_only=True)
class CladePriors:
    """These object stores priors of a single clade."""

    # The parent node denotes the most recent common ancestor (MRCA) of the
    # target clade.
    parent: UUID = field()

    # The field include the priors of all kmers shared among members of the
    # target clade.
    priors: DefaultDict[str, float] = field()


@frozen(kw_only=True)
class IngroupCladePriors:
    """These object stores priors of a single clade."""

    # The parent node denotes the most recent common ancestor (MRCA) of the
    # target clade.
    parent: UUID = field()

    # Priors of the ingroup.
    ingroup_priors: DefaultDict[str, float] = field()

    # Priors of the sister group.
    sister_group_priors: DefaultDict[str, float] = field()

    # Priors of the noise group.
    noise_group_priors: DefaultDict[str, float] = field()


@frozen(kw_only=True)
class TreePriors:
    """These object stores priors for the complete phylogenetic tree. Here nodes
    are refereed trough the linear tree clade IDs.
    """

    outgroup: CladePriors = field()
    ingroups: List[IngroupCladePriors] = field(default=[])
