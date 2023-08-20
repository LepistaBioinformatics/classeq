from copy import deepcopy
from typing import Any, Self
from uuid import UUID, uuid4

from Bio.Phylo.BaseTree import Clade, Tree


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


class ExtendedBioPythonClade(ExtendedBioPythonElement, Clade):
    # ? ------------------------------------------------------------------------
    # ? Class attributes
    # ? ------------------------------------------------------------------------

    _id: UUID
    _is_outgroup: bool

    # ? ------------------------------------------------------------------------
    # ? Life cycle hook methods
    # ? ------------------------------------------------------------------------

    def __init__(
        self,
        id: UUID | None = None,
        is_outgroup: bool | None = None,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        super().__init__(*args, **kwargs)
        self._id = id or uuid4()
        self._is_outgroup = is_outgroup or False

    # ? ------------------------------------------------------------------------
    # ? Public instance methods
    # ? ------------------------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self._id,
            "is_outgroup": self._is_outgroup,
            "branch_length": self.branch_length,
            "name": self.name,
            "confidence": self.confidence,
            "color": self.color,
            "width": self.width,
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
            is_outgroup=content.get("is_outgroup"),
            branch_length=content.get("branch_length"),
            name=content.get("name"),
            confidence=content.get("confidence"),
            color=content.get("color"),
            width=content.get("width"),
            clades=[
                cls.from_dict(clade) for clade in content.get("clades", [])
            ],
        )

    @classmethod
    def from_bio_python_clade(
        cls,
        clade: Clade,
        outgroups: list[str],
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
                    is_outgroup=clade.name in outgroups,
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
            is_outgroup=clade.name in outgroups,
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

    def set_clade_name(self, clade_id: UUID, name: str) -> None:
        for clade in [
            clade for clade in self.find_clades() if clade._id == clade_id
        ]:
            setattr(clade, "name", name)

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
    def from_bio_python_tree(cls, tree: Tree, outgroups: list[str]) -> Self:
        inner_tree = deepcopy(tree)

        return cls(
            root=ExtendedBioPythonClade.from_bio_python_clade(
                inner_tree.root,
                outgroups=outgroups,
            ),
            rooted=inner_tree.rooted,
            id=inner_tree.id,
            name=inner_tree.name,
        )
