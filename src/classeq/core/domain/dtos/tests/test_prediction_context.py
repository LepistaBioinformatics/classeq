from io import StringIO
from os import getenv
from pathlib import Path
from tempfile import NamedTemporaryFile
from unittest import TestCase

from yaml import safe_load

from classeq.core.domain.dtos.prediction_context import PredictionContext


class PredictionContextTest(TestCase):
    def setUp(self) -> None:
        INDICES_PATH = getenv("INDICES_PATH")
        ANNOTATED_PHYLOJSON_PATH = getenv("ANNOTATED_PHYLOJSON_PATH")

        if INDICES_PATH is None:
            self.skipTest(
                "INDICES_PATH is not set. Please, set it to run this test"
            )

        if ANNOTATED_PHYLOJSON_PATH is None:
            self.skipTest(
                "ANNOTATED_PHYLOJSON_PATH is not set. Please, set it to run this test"
            )

        context_content = f"""
contexts:
  - name: default
    indices: {INDICES_PATH}
    annotations: {ANNOTATED_PHYLOJSON_PATH}
    metadata:
        description: "This is a test prediction context"
"""

        content = safe_load(StringIO(context_content))
        if content is None:
            raise ValueError("Content is None")

        self.__tmp_context_file = NamedTemporaryFile(
            prefix="context", suffix=".yaml"
        )

        with open(self.__tmp_context_file.name, "w") as file:
            file.write(context_content)

    def tearDown(self) -> None:
        self.__tmp_context_file.close()

    def test_context_loading(self) -> None:
        response = PredictionContext.load_from_yaml(
            path=Path(self.__tmp_context_file.name)
        )

        self.assertFalse(response.is_left)
        self.assertTrue(response.is_right)
        self.assertIsInstance(response.value, PredictionContext)

        self.assertEqual(len(response.value.contexts), 1)
        self.assertIn("default", response.value.contexts)
        self.assertEqual(
            response.value.contexts["default"].name,
            "default",
        )

        self.assertEqual(
            response.value.contexts["default"].indices,
            Path(getenv("INDICES_PATH")),
        )

        self.assertEqual(
            response.value.contexts["default"].annotations,
            Path(getenv("ANNOTATED_PHYLOJSON_PATH")),
        )

        self.assertEqual(
            response.value.contexts["default"].metadata,
            {"description": "This is a test prediction context"},
        )
