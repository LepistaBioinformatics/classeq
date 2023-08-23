from typing import Self


class OrderedTuple(tuple[int, ...]):
    def __new__(cls, values: tuple[int, ...] | list[int] | set[int]) -> Self:
        return super(OrderedTuple, cls).__new__(cls, tuple(sorted(values)))  # type: ignore
