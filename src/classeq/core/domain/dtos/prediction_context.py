import gzip
import tarfile
from datetime import datetime
from io import BytesIO
from json import load, loads
from pathlib import Path, PosixPath
from typing import Any, Self
from uuid import UUID, uuid3, uuid4

import clean_base.exceptions as c_exc
from attrs import define, field
from clean_base.either import Either, right
from yaml import safe_load

from classeq.core.domain.dtos.biopython_wrappers import (
    ExtendedBioPythonClade,
    ExtendedBioPythonTree,
)
from classeq.core.domain.dtos.priors import TreePriors
from classeq.core.domain.dtos.reference_set import ReferenceSet
from classeq.core.domain.dtos.tests.test_reference_set import ReferenceSetTest
from classeq.settings import (
    LOGGER,
    REFERENCE_SET_OUTPUT_FILE_NAME,
    TRAIN_SOURCE_OUTPUT_FILE_NAME,
)


@define(kw_only=True)
class Context:
    name: str = field()
    indices: Path = field()
    annotations: Path = field()
    metadata: dict[str, Any] = field(default={})

    prediction_priors: tuple[ReferenceSetTest, TreePriors] = field(init=False)
    resolved_names: dict[UUID, ExtendedBioPythonClade] = field(init=False)
    id: UUID = field(init=False, factory=uuid4)

    def __attrs_post_init__(self) -> None:
        # ? --------------------------------------------------------------------
        # ? Load annotations
        # ? --------------------------------------------------------------------

        resolved_names = dict()
        if self.annotations is not None:
            with self.annotations.open("r") as f:
                annotated_tree = ExtendedBioPythonTree.from_dict(
                    content=load(f)
                )

                resolved_names = {
                    node._id: node
                    for node in annotated_tree.get_nonterminals()
                    if node.name is not None
                }

        self.resolved_names = resolved_names

        # ? --------------------------------------------------------------------
        # ? Load indices from artifact
        # ? --------------------------------------------------------------------

        if self.indices is None:
            raise ValueError("Indices path is None")

        if not self.indices.exists():
            raise ValueError(f"Indices path {self.indices} does not exist")

        if not self.indices.is_file():
            raise ValueError(f"Indices path {self.indices} is not a file")

        if (
            prediction_priors := self.__load_indices_from_classeq_artifact(
                artifact_path=self.indices
            )
        ).is_left:
            raise ValueError(prediction_priors.value.msg)

        self.prediction_priors = prediction_priors.value
        self.id = uuid3(UUID(int=0), self.name)

        LOGGER.info(
            f"Loaded context `{self.name}` with `{len(self.resolved_names)}` annotated nodes"
        )

    @staticmethod
    def __load_indices_from_classeq_artifact(
        artifact_path: Path,
    ) -> Either[c_exc.UseCaseError, tuple[ReferenceSetTest, TreePriors]]:
        try:
            with tarfile.open(artifact_path) as tf:
                # ? ------------------------------------------------------------
                # ? Load reference set
                # ? ------------------------------------------------------------

                if (
                    reference_set_content := tf.extractfile(
                        member=REFERENCE_SET_OUTPUT_FILE_NAME
                    )
                ) is None:
                    return c_exc.UseCaseError(
                        f"Reference set not found: {REFERENCE_SET_OUTPUT_FILE_NAME}",
                        logger=LOGGER,
                    )()

                with gzip.GzipFile(
                    fileobj=BytesIO(reference_set_content.read()),
                    mode="rb",
                ) as fin:
                    json_bytes = fin.read()
                    json_str = json_bytes.decode("utf-8")

                if (
                    reference_set_either := ReferenceSet.from_dict(
                        content=loads(json_str)
                    )
                ).is_left:
                    return reference_set_either

                # ? ------------------------------------------------------------
                # ? Load indices
                # ? ------------------------------------------------------------

                if (
                    indices_content := tf.extractfile(
                        member=TRAIN_SOURCE_OUTPUT_FILE_NAME
                    )
                ) is None:
                    return c_exc.UseCaseError(
                        f"Indices not found: {TRAIN_SOURCE_OUTPUT_FILE_NAME}",
                        logger=LOGGER,
                    )()

                with gzip.GzipFile(
                    fileobj=BytesIO(indices_content.read()),
                    mode="rb",
                ) as fin:
                    json_bytes = fin.read()
                    json_str = json_bytes.decode("utf-8")

                if (
                    tree_priors_either := TreePriors.from_dict(
                        content=loads(json_str)
                    )
                ).is_left:
                    return tree_priors_either

            return right((reference_set_either.value, tree_priors_either.value))

        except Exception as exc:
            return c_exc.UseCaseError(exc, logger=LOGGER)()


@define(kw_only=True)
class PredictionContext:
    last_reload: datetime = field(default=datetime.now())
    contexts: dict[UUID, Context] = field(default={})

    def update_last_reload(self) -> None:
        self.last_reload = datetime.now()

    def get_context_by_name(
        self,
        identifier: UUID,
    ) -> Either[c_exc.UseCaseError, Context | None]:
        return self.contexts.get(identifier)

    @classmethod
    def load_from_yaml(cls, path: Path) -> Either[c_exc.UseCaseError, Self]:
        try:
            if not isinstance(path, PosixPath):
                return c_exc.UseCaseError(
                    f"Path must be an instance of Path, got {type(path)}",
                    logger=LOGGER,
                )()

            if not path.exists():
                return c_exc.UseCaseError(
                    f"Path {path} does not exist",
                    logger=LOGGER,
                )()

            with path.open("r") as file:
                content = safe_load(file)

            if content is None:
                return c_exc.UseCaseError(
                    "Content is None",
                    logger=LOGGER,
                )()

            if (contexts := content.get("contexts")) is None:
                return c_exc.UseCaseError(
                    "Invalid content, missing 'contexts' key.",
                    logger=LOGGER,
                )()

            if not isinstance(contexts, list):
                return c_exc.UseCaseError(
                    "Invalid content, 'contexts' key must be a list.",
                    logger=LOGGER,
                )()

            if not contexts:
                return c_exc.UseCaseError(
                    "Invalid content, 'contexts' key must not be empty.",
                    logger=LOGGER,
                )()

            if len(set([c.get("name") for c in contexts])) != len(contexts):
                return c_exc.UseCaseError(
                    "Invalid content, 'contexts' key must not contain duplicated names.",
                    logger=LOGGER,
                )()

            output_contexts: dict[UUID, Context] = dict()
            for context in contexts:
                if (name := context.get("name")) is None:
                    return c_exc.UseCaseError(
                        "Invalid content, missing 'name' key.",
                        logger=LOGGER,
                    )()

                if not isinstance(name, str):
                    return c_exc.UseCaseError(
                        "Invalid content, 'name' key must be a string.",
                        logger=LOGGER,
                    )()

                if (indices := context.get("indices")) is None:
                    return c_exc.UseCaseError(
                        "Invalid content, missing 'indices' key.",
                        logger=LOGGER,
                    )()

                if not isinstance(indices, str):
                    return c_exc.UseCaseError(
                        "Invalid content, 'indices' key must be a string.",
                        logger=LOGGER,
                    )()

                if (annotations := context.get("annotations")) is None:
                    return c_exc.UseCaseError(
                        "Invalid content, missing 'annotations' key.",
                        logger=LOGGER,
                    )()

                if not isinstance(annotations, str):
                    return c_exc.UseCaseError(
                        "Invalid content, 'annotations' key must be a string.",
                        logger=LOGGER,
                    )()

                if (metadata := context.get("metadata", {})) is not None:
                    if not isinstance(metadata, dict):
                        return c_exc.UseCaseError(
                            "Invalid content, 'metadata' key must be a dict.",
                            logger=LOGGER,
                        )()

                new_context = Context(
                    name=name,
                    indices=Path(indices),
                    annotations=Path(annotations),
                    metadata=metadata,
                )

                output_contexts.update({new_context.id: new_context})

            instance = cls(contexts=output_contexts)
            instance.update_last_reload()
            return right(instance)

        except Exception as exc:
            return c_exc.UseCaseError(exc, logger=LOGGER)()
