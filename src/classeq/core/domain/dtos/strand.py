from enum import Enum
from typing import Any


class StrandEnum(Enum):
    BOTH = "both"
    PLUS = "plus"
    MINUS = "minus"

    @classmethod
    def _missing_(cls, _: object) -> Any:
        return cls.BOTH
