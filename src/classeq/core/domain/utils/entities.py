from __future__ import annotations

from abc import ABCMeta, abstractmethod
from typing import Any, NamedTuple, NoReturn

from .either import Either
from .exceptions import (
    CreationError,
    DeletionError,
    ExecutionError,
    FetchingError,
    UpdatingError,
)

# ------------------------------------------------------------------------------
# DEFAULT RESPONSE TYPES
# ------------------------------------------------------------------------------


class CreateResponse(NamedTuple):
    """A interface for creation return methods."""

    created: bool
    instance: Any


class GetOrCreateResponse(NamedTuple):
    """A interface for get_or_creation return methods."""

    created: bool
    instance: Any


class FetchResponse(NamedTuple):
    """A interface for fetch return methods."""

    fetched: bool
    instance: NoReturn | Any


class FetchManyResponse(NamedTuple):
    """A interface for fetch return methods."""

    fetched: bool
    instance: NoReturn | list[Any]


class UpdateResponse(NamedTuple):
    """A interface for update return methods."""

    updated: bool
    instance: Any


class UpdateManyResponse(NamedTuple):
    """A interface for update return methods."""

    updated: bool
    records: int


class DeleteResponse(NamedTuple):
    """A interface for delete return methods."""

    deleted: bool
    instance: Any


class ExecuteResponse(NamedTuple):
    """A interface for execution return methods."""

    executed: bool
    instance: Any


class Fetching(metaclass=ABCMeta):
    """The registration interface."""

    # --------------------------------------------------------------------------
    # CLASS ATTRIBUTES
    # --------------------------------------------------------------------------

    tclass: Any

    # --------------------------------------------------------------------------
    # MAGIC METHODS
    # --------------------------------------------------------------------------

    def __init_subclass__(cls, tclass: Any) -> None:
        cls.tclass = tclass

    # --------------------------------------------------------------------------
    # ABSTRACT METHODS
    # --------------------------------------------------------------------------

    @abstractmethod
    def show(
        self,
        term: str | None = None,
        page: int = 1,
        size: int = 10,
        **kwargs: Any,
    ) -> Either[FetchingError, FetchManyResponse]:
        raise NotImplementedError

    @abstractmethod
    def get(
        self,
        pk: int | None = None,
        other: str | None = None,
        **kwargs: Any,
    ) -> Either[FetchingError, FetchResponse]:
        raise NotImplementedError


class Registration(metaclass=ABCMeta):
    """The registration interface."""

    # --------------------------------------------------------------------------
    # CLASS ATTRIBUTES
    # --------------------------------------------------------------------------

    tclass: Any

    # --------------------------------------------------------------------------
    # MAGIC METHODS
    # --------------------------------------------------------------------------

    def __init_subclass__(cls, tclass: Any) -> None:
        cls.tclass = tclass

    # --------------------------------------------------------------------------
    # ABSTRACT METHODS
    # --------------------------------------------------------------------------

    @abstractmethod
    def create(
        self, registration_type: Registration.tclass, **kwargs: Any
    ) -> Either[CreationError, CreateResponse]:
        raise NotImplementedError

    @abstractmethod
    def get_or_create(
        self, registration_type: Registration.tclass, **kwargs: Any
    ) -> Either[CreationError, GetOrCreateResponse]:
        raise NotImplementedError


class Updating(metaclass=ABCMeta):
    """The updating interface."""

    # --------------------------------------------------------------------------
    # CLASS ATTRIBUTES
    # --------------------------------------------------------------------------

    tclass: Any

    # --------------------------------------------------------------------------
    # MAGIC METHODS
    # --------------------------------------------------------------------------

    def __init_subclass__(cls, tclass: Any) -> None:
        cls.tclass = tclass

    # --------------------------------------------------------------------------
    # ABSTRACT METHODS
    # --------------------------------------------------------------------------

    @abstractmethod
    def update(
        self, updating_type: Updating.tclass, **kwargs: Any
    ) -> Either[UpdatingError, UpdateResponse]:
        raise NotImplementedError


class Deletion(metaclass=ABCMeta):
    """The deletion interface."""

    # --------------------------------------------------------------------------
    # CLASS ATTRIBUTES
    # --------------------------------------------------------------------------

    tclass: Any

    # --------------------------------------------------------------------------
    # MAGIC METHODS
    # --------------------------------------------------------------------------

    def __init_subclass__(cls, tclass: Any) -> None:
        cls.tclass = tclass

    # --------------------------------------------------------------------------
    # ABSTRACT METHODS
    # --------------------------------------------------------------------------

    @abstractmethod
    def delete(
        self, deletion_type: Deletion.tclass, **kwargs: Any
    ) -> Either[DeletionError, DeleteResponse]:
        raise NotImplementedError


class ExecuteStep(metaclass=ABCMeta):
    """The execution interface."""

    # --------------------------------------------------------------------------
    # CLASS ATTRIBUTES
    # --------------------------------------------------------------------------

    tclass: Any

    # --------------------------------------------------------------------------
    # MAGIC METHODS
    # --------------------------------------------------------------------------

    def __init_subclass__(cls, tclass: Any) -> None:
        cls.tclass = tclass

    # --------------------------------------------------------------------------
    # ABSTRACT METHODS
    # --------------------------------------------------------------------------

    @abstractmethod
    def run(
        self, execution_type: ExecuteStep.tclass, **kwargs: Any
    ) -> Either[ExecutionError, ExecuteResponse]:
        raise NotImplementedError


# ------------------------------------------------------------------------------
# SETUP DEFAULT EXPORTS
# ------------------------------------------------------------------------------


__all__ = [
    "CreateResponse",
    "GetOrCreateResponse",
    "FetchResponse",
    "FetchManyResponse",
    "UpdateResponse",
    "UpdateManyResponse",
    "DeleteResponse",
    "ExecuteResponse",
    "Fetching",
    "Registration",
    "Updating",
    "Deletion",
    "ExecuteStep",
]
