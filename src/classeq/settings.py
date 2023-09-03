from logging import (
    CRITICAL,
    DEBUG,
    ERROR,
    FATAL,
    INFO,
    NOTSET,
    WARN,
    WARNING,
    FileHandler,
    Formatter,
    StreamHandler,
    getLogger,
)
from os import getenv

# ? ----------------------------------------------------------------------------
# ? Build logger configurations
# ? ----------------------------------------------------------------------------


LOGGING_LEVEL = INFO


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


LOGGER = getLogger("CLASSEQ")
LOGGER.setLevel(DEBUG)


formatter = Formatter(
    "%(levelname)s [ %(asctime)s ] %(message)s",
    datefmt="%Y-%d-%m %H:%M:%S",
)


file_handler = FileHandler(getenv("LOGGING_FILE", "classeq.log"))
file_handler.setLevel(DEBUG)
file_handler.setFormatter(formatter)
LOGGER.addHandler(file_handler)


stream_handler = StreamHandler()
stream_handler.setLevel(LOGGING_LEVEL)
stream_handler.setFormatter(formatter)
LOGGER.addHandler(stream_handler)


# ? ----------------------------------------------------------------------------
# ? DNA bases used o kmer count
# ? ----------------------------------------------------------------------------


BASES = ["A", "C", "T", "G"]


DEFAULT_KMER_SIZE = 12


# ? ----------------------------------------------------------------------------
# ? Minimum clade size
# ? ----------------------------------------------------------------------------


MINIMUM_CLADE_SIZE = 2


# ? ----------------------------------------------------------------------------
# ? Prediction related params
# ? ----------------------------------------------------------------------------


MINIMUM_INGROUP_QUERY_KMERS_MATCH = 3


MINIMUM_INGROUP_SISTER_MATCH_KMERS_DIFFERENCE = 2


# ? ----------------------------------------------------------------------------
# ? Temp files names
# ? ----------------------------------------------------------------------------


TEMP_INPUT_FILE_SUFFIX = "sanitized"


DEFAULT_CLASSEQ_OUTPUT_FILE_NAME = "classeq.tar"


REFERENCE_SET_OUTPUT_FILE_NAME = "reference-set.json.gz"


TRAIN_SOURCE_OUTPUT_FILE_NAME = "train-source.json.gz"
