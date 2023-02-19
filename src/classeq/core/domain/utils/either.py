from typing import Any, Generic, TypeVar, Union

_L = TypeVar("_L")
_R = TypeVar("_R")


class __Base:
    """A superclass of `Left` and `Right`."""

    _is_left: bool = True
    _is_right: bool = False
    _value: Any

    @property
    def is_left(self) -> bool:
        return self._is_left

    @property
    def is_right(self) -> bool:
        return self._is_right


class Left(__Base, Generic[_L]):
    """A class to transport undesired responses."""

    _is_left: bool = True
    _is_right: bool = False

    def __init__(self, value: _L) -> None:
        self._value = value

    @property
    def value(self) -> _L:
        return self._value


class Right(__Base, Generic[_R]):
    """A class to transport right responses."""

    _is_left: bool = False
    _is_right: bool = True

    def __init__(self, value: _R) -> None:
        self._value = value

    @property
    def value(self) -> _R:
        return self._value


Either = Union[Left[_L], Right[_R]]


def left(left: _L) -> Left[_L]:
    """Create a Left instance.

    Args:
        left (_L): The value as a left content.

    Returns:
        Left: A Left instance.
    """

    return Left(left)


def right(right: _R) -> Right[_R]:
    """Create a Right instance.

    Args:
        right (_R): The value as a right content.

    Returns:
        Right: A Right instance.
    """

    return Right(right)


# ------------------------------------------------------------------------------
# SETUP DEFAULT EXPORTS
# ------------------------------------------------------------------------------


__all__ = ["Left", "Right", "Either", "left", "right"]
