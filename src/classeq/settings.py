from logging import (
    CRITICAL,
    DEBUG,
    ERROR,
    FATAL,
    INFO,
    NOTSET,
    WARN,
    WARNING,
    getLogger,
)
from os import getenv

from classeq.core.domain.utils.validations import should_be_int


# ? ----------------------------------------------------------------------------
# ? Build logger configurations
# ? ----------------------------------------------------------------------------


LOGGING_LEVEL: int = DEBUG


ENV_LOGGING_LEVEL = getenv("LOGGING_LEVEL")


if ENV_LOGGING_LEVEL is not None:
    ENV_LOGGING_LEVEL = ENV_LOGGING_LEVEL.upper()

    if should_be_int(ENV_LOGGING_LEVEL) is False:
        raise ValueError(
            f"Invalid `ENV_LOGGING_LEVEL` type: {ENV_LOGGING_LEVEL}"
        )

    ENV_LOGGING_LEVEL = int(ENV_LOGGING_LEVEL)  # type: ignore

    if ENV_LOGGING_LEVEL not in [
        CRITICAL,
        DEBUG,
        ERROR,
        FATAL,
        INFO,
        NOTSET,
        WARN,
        WARNING,
    ]:
        raise ValueError(
            f"Invalid `ENV_LOGGING_LEVEL` value from env: {ENV_LOGGING_LEVEL}"
        )

    LOGGING_LEVEL = ENV_LOGGING_LEVEL  # type: ignore


# ? ----------------------------------------------------------------------------
# ? Initialize global logger
# ? ----------------------------------------------------------------------------


LOGGER = getLogger("CLASSEQ")


LOGGER.setLevel(LOGGING_LEVEL)


# ? ----------------------------------------------------------------------------
# ? DNA bases used o kmer count
# ? ----------------------------------------------------------------------------


BASES = ["A", "C", "T", "G"]


DEFAULT_KMER_SIZE = 8
