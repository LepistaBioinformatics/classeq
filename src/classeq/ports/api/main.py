from contextlib import asynccontextmanager
from os import getenv
from pathlib import Path
from typing import Any, Callable
from uuid import UUID

from Bio.SeqRecord import SeqRecord
from fastapi import Depends, FastAPI, HTTPException, status
from pydantic import BaseModel
from starlette.middleware.cors import CORSMiddleware

from classeq.core.domain.dtos.prediction_context import PredictionContext
from classeq.core.use_cases.build_prediction_context import (
    build_prediction_context,
)
from classeq.core.use_cases.build_prediction_context.errors import (
    FailLoadingContext,
    InvalidReferenceTreeFromContext,
    UnexistingContext,
    UnexpectedErrorOnExecuteAdherenceTest,
)
from classeq.core.use_cases.predict_taxonomies_recursively._perform_single_sequence_phylogenetic_adherence_test import (
    PredictionResult,
)
from classeq.core.use_cases.predict_taxonomies_recursively._perform_adherence_test_of_child_clades import (
    CladeAdherenceResultStatus,
)
from classeq.settings import LOGGER

# ? ----------------------------------------------------------------------------
# ? Context configuration
# ? ----------------------------------------------------------------------------


__PREDICTION_CONTEXT: tuple[
    PredictionContext, Callable[..., PredictionResult]
] | None = None


# ? ----------------------------------------------------------------------------
# ? API configuration
# ? ----------------------------------------------------------------------------


@asynccontextmanager
async def lifespan(_: FastAPI) -> Any:
    global __PREDICTION_CONTEXT

    CONTEXT_ENV_VAR = getenv("CONTEXT_MODELS", None)

    if CONTEXT_ENV_VAR is None:
        raise ValueError("Missing `CONTEXT_MODELS` environment variable")

    if (
        context_either := build_prediction_context(
            context_path=Path(CONTEXT_ENV_VAR)
        )
    ).is_left:
        raise ValueError(context_either.value)

    __PREDICTION_CONTEXT = context_either.value

    yield


__API_BASE_URL = "/api"


app = FastAPI(
    title="Classeq API",
    openapi_url=f"{__API_BASE_URL}/openapi.json",
    docs_url=f"{__API_BASE_URL}/docs",
    lifespan=lifespan,
)


app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["GET", "POST"],
    allow_headers=["*"],
)


# ? ----------------------------------------------------------------------------
# ? Auxiliary functions
# ? ----------------------------------------------------------------------------


async def __get_predictor() -> Callable[..., PredictionResult]:
    if __PREDICTION_CONTEXT is None:
        raise HTTPException(
            detail="No models loaded",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )

    return __PREDICTION_CONTEXT[1]


# ? ----------------------------------------------------------------------------
# ? API Models
# ? ----------------------------------------------------------------------------


class ModelResponse(BaseModel):
    id: UUID
    name: str
    metadata: dict[str, Any]


class Query(BaseModel):
    header: str
    sequence: str


class PredictionStepResponse(BaseModel):
    depth: int

    # TODO:
    #
    # Do implement a pydantiv representation of the original
    # `CladeAdherenceResult` class as a response model.
    result: Any


class PredictionResultResponse(BaseModel):
    status: CladeAdherenceResultStatus
    path: list[PredictionStepResponse]


class PredictionResponse(BaseModel):
    model: UUID
    query: str
    prediction: PredictionResultResponse


# ? ----------------------------------------------------------------------------
# ? API Endpoints
# ? ----------------------------------------------------------------------------


@app.get("/api/models/")
async def get_models() -> list[ModelResponse]:
    if __PREDICTION_CONTEXT is None:
        raise HTTPException(
            detail="No models loaded",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )

    return [
        ModelResponse(
            **{
                "id": model,
                "name": value.name,
                "metadata": value.metadata,
            }
        )
        for model, value in __PREDICTION_CONTEXT[0].contexts.items()
    ]


@app.post("/api/predict/{id}")
async def predict(
    id: UUID,
    query: Query,
    predictor: Callable[..., PredictionResult] = Depends(__get_predictor),
) -> PredictionResponse:
    try:
        prediction: PredictionResult = predictor(
            model_id=id,
            seq=SeqRecord(
                id=query.header,
                seq=query.sequence,
            ),
        )

    except UnexistingContext as e:
        LOGGER.exception(e)
        raise HTTPException(
            detail=str(e),
            status_code=status.HTTP_400_BAD_REQUEST,
        )

    except (
        FailLoadingContext,
        InvalidReferenceTreeFromContext,
        UnexpectedErrorOnExecuteAdherenceTest,
    ) as e:
        LOGGER.exception(e)
        raise HTTPException(
            detail=str(e),
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )

    except Exception as e:
        LOGGER.exception(e)
        raise HTTPException(
            detail="Unexpected error",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        )

    return PredictionResponse(
        model=id,
        query=query.header,
        prediction=prediction.to_dict(omit_children=True),
    )
