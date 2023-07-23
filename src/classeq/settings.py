from logging import (
    CRITICAL,
    DEBUG,
    ERROR,
    FATAL,
    INFO,
    NOTSET,
    WARN,
    WARNING,
    basicConfig,
    getLogger,
)
from os import getenv

# ? ----------------------------------------------------------------------------
# ? Build logger configurations
# ? ----------------------------------------------------------------------------


LOGGING_LEVEL = DEBUG


ENV_LOGGING_LEVEL = getenv("LOGGING_LEVEL")


if ENV_LOGGING_LEVEL is not None:
    ENV_LOGGING_LEVEL = eval(ENV_LOGGING_LEVEL.upper())

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


basicConfig(
    level=DEBUG,
    format="%(levelname)s\t[ %(asctime)s ]\t%(message)s",
)


LOGGER = getLogger("CLASSEQ")


LOGGER.setLevel(LOGGING_LEVEL)


# ? ----------------------------------------------------------------------------
# ? DNA bases used o kmer count
# ? ----------------------------------------------------------------------------


BASES = ["A", "C", "T", "G"]


DEFAULT_KMER_SIZE = 32


# ? ----------------------------------------------------------------------------
# ? Temp files names
# ? ----------------------------------------------------------------------------


TEMP_INPUT_FILE_SUFFIX = "sanitized"
