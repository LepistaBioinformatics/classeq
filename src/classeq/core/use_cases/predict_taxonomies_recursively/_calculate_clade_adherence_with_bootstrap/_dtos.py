from enum import Enum

from attrs import field, frozen


class AdherenceStatus(Enum):
    SUCCESS = "success"
    NOT_ENOUGH_PRIORS = "not_enough_priors"
    UNDEFINED = "undefined"


@frozen(kw_only=True)
class AdherenceResult:
    used_kmers: int = field()
    query_kmers_size: int = field()
    subject_kmers_size: int = field()
    status: AdherenceStatus = field(default=AdherenceStatus.UNDEFINED)
    joint_probability: float = field(default=-999)

    def __str__(self) -> str:
        return (
            "AdherenceResult( "
            + f"used={self.used_kmers}, "
            + f"q={self.query_kmers_size}, "
            + f"s={self.subject_kmers_size}, "
            + f"status={self.status.name}, "
            + "joint_probability=%0.10f )" % self.joint_probability
        )
