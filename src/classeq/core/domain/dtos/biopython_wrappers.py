from copy import deepcopy
from enum import Enum
from typing import Any, Self
from uuid import UUID, uuid3

from Bio.Phylo.BaseTree import Clade, Tree

from classeq.settings import CLASSEQ_NAMESPACE_DNS


class ExtendedBioPythonElement:
    # ? ------------------------------------------------------------------------
    # ? Life cycle hook methods
    # ? ------------------------------------------------------------------------

    def __repr__(self) -> str:
        """Show this object's constructor with its primitive arguments."""

        def pair_as_kwarg_string(key: str, val: Any) -> str:
            if isinstance(val, str):
                val = val[:57] + "..." if len(val) > 60 else val
                return f"{key}='{val}'"
            return f"{key}={val}"

        return "%s(%s)" % (
            self.__class__.__name__,
            ", ".join(
                pair_as_kwarg_string(key, val)
                for key, val in sorted(self.__dict__.items())
                if val is not None
                and type(val) in (str, int, float, bool, str, UUID)
            ),
        )

    __str__ = __repr__  # type: ignore


class MajorTaxonomicRanks(Enum):
    KINGDOM = "kingdom"
    PHYLUM = "phylum"
    CLASS = "class"
    ORDER = "order"
    FAMILY = "family"
    GENUS = "genus"
    SPECIES = "species"


class MinorTaxonomicRanks(Enum):
    SUB_KINGDOM = "subkingdom"
    SUB_PHYLUM = "subphylum"
    INFRA_CLASS = "infraclass"
    SUB_CLASS = "subclass"
    SUB_ORDER = "suborder"
    INFRA_ORDER = "infraorder"
    SUB_FAMILY = "subfamily"
    SUB_GENUS = "subgenus"
    SUB_SPECIES = "subspecies"


class ExtraTaxonomicRanks(Enum):
    NO_RANK = "no-rank"
    BIOTYPE = "biotype"
    CLADE = "clade"
    COHORT = "cohort"
    SUB_COHORT = "subcohort"
    FORMA = "forma"
    FORMA_SPECIALIS = "forma-specialis"
    GENOTYPE = "genotype"
    ISOLATE = "isolate"
    MORPH = "morph"
    PARV_ORDER = "parvorder"
    PATHO_GROUP = "pathogroup"
    SECTION = "section"
    SERIES = "series"
    SERO_GROUP = "serogroup"
    SEROTYPE = "serotype"
    SPECIES_GROUP = "species-group"
    SPECIES_SUBGROUP = "species-subgroup"
    STRAIN = "strain"
    SUB_SECTION = "subsection"
    SUB_TRIBE = "subtribe"
    SUB_VARIETY = "subvariety"
    SUPER_CLASS = "superclass"
    SUPER_FAMILY = "superfamily"
    SUPER_KINGDOM = "superkingdom"
    SUPER_ORDER = "superorder"
    SUPER_PHYLUM = "superphylum"
    TRIBE = "tribe"
    VARIETAS = "varietas"


def try_to_reach_rank_enum(
    rank: str,
) -> MajorTaxonomicRanks | MinorTaxonomicRanks | ExtraTaxonomicRanks | None:
    try:
        return MajorTaxonomicRanks(rank)
    except ValueError:
        pass

    try:
        return MinorTaxonomicRanks(rank)
    except ValueError:
        pass

    try:
        return ExtraTaxonomicRanks(rank)
    except ValueError:
        pass

    return None


class ExtendedBioPythonClade(ExtendedBioPythonElement, Clade):
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    _id: UUID
    _taxid: int | None = None
    _related_rank: MajorTaxonomicRanks | MinorTaxonomicRanks | ExtraTaxonomicRanks | None = (
        None
    )

    # ? ------------------------------------------------------------------------
    # ? Life cycle hook methods
    # ? ------------------------------------------------------------------------

    def __init__(
        self,
        id: UUID | None = None,
        taxid: int | None = None,
        related_rank: MajorTaxonomicRanks
        | MinorTaxonomicRanks
        | ExtraTaxonomicRanks
        | None = None,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        super().__init__(*args, **kwargs)

        self._id = id or uuid3(
            CLASSEQ_NAMESPACE_DNS,
            (
                (self.name or "")
                + "".join(
                    [
                        "".join([j.name for j in i.get_terminals()])
                        for i in self.clades
                    ]
                    or []
                )
            ),
        )

        self._taxid = taxid or None
        self._related_rank = related_rank or None

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self._id,
            "branch_length": self.branch_length,
            "name": self.name,
            "confidence": self.confidence,
            "color": self.color,
            "width": self.width,
            "taxid": self._taxid,
            "related_rank": (
                self._related_rank.value
                if self._related_rank is not None
                else None
            ),
            "clades": [clade.to_dict() for clade in self.clades],
        }

    # ? ------------------------------------------------------------------------
    # ? Public class methods
    # ? ------------------------------------------------------------------------

    @classmethod
    def from_dict(
        cls,
        content: dict[str, Any],
    ) -> Self:
        return cls(
            id=UUID(content.get("id")),
            branch_length=content.get("branch_length"),
            name=content.get("name"),
            confidence=content.get("confidence"),
            color=content.get("color"),
            width=content.get("width"),
            taxid=content.get("taxid"),
            related_rank=(
                try_to_reach_rank_enum(related_rank)
                if (related_rank := content.get("related_rank")) is not None
                else None
            ),
            clades=[
                cls.from_dict(clade) for clade in content.get("clades", [])
            ],
        )

    @classmethod
    def from_bio_python_clade(
        cls,
        clade: Clade,
    ) -> Self:
        """Recursively copy a BioPython Clade object and its children.

        Args:
            clade (Clade): A BioPython Clade object.

        Returns:
            Self: A copy of the input Clade object as `ExtendedBioPythonClade`.

        """

        def recursive_parse_clades(clades: list[Clade]) -> list[Clade]:
            clades_copy: list[Clade] = []
            for clade in clades:
                copied_clade = ExtendedBioPythonClade(
                    branch_length=clade.branch_length,
                    name=clade.name,
                    clades=clade.clades,
                    confidence=clade.confidence,
                    color=clade.color,
                    width=clade.width,
                )

                copied_clade.clades = recursive_parse_clades(
                    copied_clade.clades
                )

                clades_copy.append(copied_clade)

            return clades_copy

        clade.clades = recursive_parse_clades(clade.clades)

        return cls(
            branch_length=clade.branch_length,
            name=clade.name,
            clades=clade.clades,
            confidence=clade.confidence,
            color=clade.color,
            width=clade.width,
        )


class ExtendedBioPythonTree(Tree):
    # ? ------------------------------------------------------------------------
    # ? Life cycle hook methods
    # ? ------------------------------------------------------------------------

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return {
            "root": self.root.to_dict(),
            "rooted": self.rooted,
            "id": self.id,
            "name": self.name,
        }

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def set_clade_name(
        self,
        id: UUID,
        name: str,
    ) -> None:
        for clade in [clade for clade in self.find_clades() if clade._id == id]:
            setattr(clade, "name", name)

    def set_clade_taxid(
        self,
        id: UUID,
        taxid: int,
    ) -> None:
        for clade in [clade for clade in self.find_clades() if clade._id == id]:
            setattr(clade, "_taxid", taxid)

    def set_clade_rank(
        self,
        id: UUID,
        rank: MajorTaxonomicRanks
        | MinorTaxonomicRanks
        | ExtraTaxonomicRanks
        | None = None,
    ) -> None:
        for clade in [clade for clade in self.find_clades() if clade._id == id]:
            setattr(clade, "_related_rank", rank)

    def find_clade_by_id(
        self,
        id: UUID,
    ) -> None:
        try:
            return next(
                clade for clade in self.find_clades() if clade._id == id
            )
        except StopIteration:
            return None

    # ? ------------------------------------------------------------------------
    # ? Public class methods
    # ? ------------------------------------------------------------------------

    @classmethod
    def from_dict(
        cls,
        content: dict[str, Any],
    ) -> Self:
        return cls(
            root=ExtendedBioPythonClade.from_dict(content.get("root")),  # type: ignore
            rooted=content.get("rooted"),
            id=content.get("id"),
            name=content.get("name"),
        )

    @classmethod
    def from_bio_python_tree(
        cls,
        tree: Tree,
    ) -> Self:
        inner_tree = deepcopy(tree)

        return cls(
            root=ExtendedBioPythonClade.from_bio_python_clade(
                clade=inner_tree.root
            ),
            rooted=inner_tree.rooted,
            id=inner_tree.id,
            name=inner_tree.name,
        )
