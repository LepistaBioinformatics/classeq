#!/usr/bin/env python3
#
#
# This is an integration test to predict the identity of a set of sequences
# using the classeq API.
#
# Before start the test you need to start the classeq API server. You can do it
# by running the following commands from the project root directory:
#
# ```bash
#
# export CONTEXT_MODELS="docs/manuscript/assays/bacillus-subtilis-group/supplementary-material/05-config.yaml"
#
# uvicorn classeq.ports.api.main:app \
#     --reload \
#     --host 0.0.0.0 \
#     --port 8080
#
# ```
#
# Then, you need to replace the `QUERY_FILE_PATH` environment variable with the
# path to the FASTA file containing the sequences to be predicted. Additionally,
# you can also replace the `MODEL_ID` environment variable with the ID of the
# model to be used. By default the model ID a UUID3 generated on the  moment
# which the model configuration file is loaded concatenating the name of the
# model with a simple UUID of zero. This is done to ensure that the model ID is
# always the same of the default presented in line 41 of this file.
#
# Case the default MODEL ID was not compatible with the model you want to use,
# you can discover the ID of the model simply by running the classeq API server
# and accessing the swagger documentation at http://localhost:8080/api/docs.
# Then, you can copy the model ID from the swagger documentation and paste it in
# the `MODEL_ID`.
#

from asyncio import gather
from os import getenv
from pathlib import Path
from typing import Any, Awaitable, Generator

from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord
from aiohttp import ClientSession, ClientTimeout


MODEL_ID = getenv("MODEL_ID", "61a72fcc-a485-3b2f-b3fa-7bedb24cc538")
API_URL = f"http://localhost:8080/api/predict/{MODEL_ID}"
QUERY_FILE_PATH = getenv("QUERY_FILE_PATH")


if not QUERY_FILE_PATH:
    raise ValueError("`QUERY_FILE_PATH` environment variable is not defined.")


def __load_query_sequences(
    query_file: Path,
) -> Generator[SeqRecord, None, None]:
    with open(query_file, "r") as handle:
        for record in parse(handle, "fasta"):
            yield record


async def __fetch_identity(
    seq: SeqRecord, session: ClientSession
) -> dict[str, Any]:
    async with session.post(
        API_URL,
        ssl=False,
        json=dict(
            header=seq.id,
            sequence=str(seq.seq),
        ),
    ) as response:
        status = response.status
        data = await response.json()

        return {"status": status, "response": data}


async def __fetch_identities() -> Awaitable[list[dict[str, Any]]]:
    async with ClientSession(
        timeout=ClientTimeout(total=60 * 60 * 24),
    ) as session:
        tasks = [
            __fetch_identity(
                seq=record,
                session=session,
            )
            for record in __load_query_sequences(
                query_file=Path(QUERY_FILE_PATH),
            )
        ]

        return await gather(*tasks)


if __name__ == "__main__":
    from asyncio import run

    run(__fetch_identities())
