import sys
from json import dump, load
from pathlib import Path
from typing import Any, Generator, Literal
from uuid import UUID

import clean_base.exceptions as c_exc
from clean_base.either import Either, right
from PySide2 import QtCore, QtGui
from PySide2.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QFileDialog,
    QInputDialog,
    QLineEdit,
    QMainWindow,
    QMenu,
    QPushButton,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from classeq.core.domain.dtos.biopython_wrappers import (
    ExtendedBioPythonClade,
    ExtendedBioPythonTree,
)
from classeq.core.domain.dtos.tree_source_format import TreeSourceFormatEnum
from classeq.ports.app.settings import STYLE_SHEET_DARK
from classeq.settings import LOGGER


class TreeEditor(QMainWindow):
    # ? ------------------------------------------------------------------------
    # ? CLASS ATTRIBUTES
    # ? ------------------------------------------------------------------------

    __tree: ExtendedBioPythonTree | None = None
    __tree_file_path: Path

    # ? ------------------------------------------------------------------------
    # ? LIFE CYCLE HOOKS
    # ? ------------------------------------------------------------------------

    def __init__(
        self,
        tree_file_path: Path,
        parent: Any = None,
    ) -> None:
        super(TreeEditor, self).__init__(parent)

        self.setWindowTitle("Classeq Tree Editor")
        self.setGeometry(QtCore.QRect(200, 200, 800, 600))
        self.setMinimumHeight(150)

        self.__main = QWidget(self)
        self.setCentralWidget(self.__main)

        self.__add_widgets()
        self.__build_main_menu(self)

        self.__tree_file_path = tree_file_path

        if (
            load_tree_either := self.__load_tree_from_phylojson_file(
                tree_file_path
            )
        ).is_left:
            raise Exception(load_tree_either.value.msg)

    # ? ------------------------------------------------------------------------
    # ? PRIVATE STATIC METHODS
    # ? ------------------------------------------------------------------------

    def __add_widgets(self) -> None:
        # ? --------------------------------------------------------------------
        # ? Main layout
        #
        # Here the phylogeny is presented.
        #
        # ? --------------------------------------------------------------------
        primary_hbox = QVBoxLayout(self.__main)

        # ? --------------------------------------------------------------------
        # ? Filter box
        # ? --------------------------------------------------------------------

        self.__search_field = QLineEdit()
        self.__search_field.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.__search_field.setObjectName("search_field")
        self.__search_field.setStyleSheet("font-size: 16px; height: 50px;")
        self.__search_field.setPlaceholderText("Search for a node")

        self.__search_field.textChanged.connect(
            lambda text: self.__on_text_changed(text)
        )

        primary_hbox.addWidget(self.__search_field)

        # ? --------------------------------------------------------------------
        # ? Internal nodes box
        # ? --------------------------------------------------------------------
        inner_tree_items = self.__nodes_tree_widget = QTreeWidget()
        inner_tree_items.setHeaderItem(
            QTreeWidgetItem(
                [
                    "Name",
                    "Support",
                    "Outgroup",
                    "Terminal",
                    "Identifier",
                    "",
                ]
            )
        )

        inner_tree_items.setColumnCount(6)
        inner_tree_items.setColumnWidth(0, 700)
        inner_tree_items.setColumnWidth(1, 80)
        inner_tree_items.setColumnWidth(2, 80)
        inner_tree_items.setColumnWidth(3, 80)
        inner_tree_items.setColumnWidth(4, 300)
        inner_tree_items.setColumnWidth(5, 50)
        inner_tree_items.setSelectionMode(QAbstractItemView.ExtendedSelection)
        inner_tree_items.setSortingEnabled(True)
        primary_hbox.addWidget(inner_tree_items)

        return

    def __on_text_changed(self, text: str) -> None:
        self.__nodes_tree_widget.clearSelection()

        if text == "" or text is None or text.__len__() < 3:
            return

        name_items: list[QTreeWidgetItem] = []
        names_field_index = 0
        ids_field_index = 4

        def recursive_find_in_children(
            items: list[QTreeWidgetItem],
            text: str,
            field_index: int,
        ) -> Generator[list[QTreeWidgetItem], None, None]:
            """Recursively find the nodes that match the search criteria.

            Args:
                items (list[QTreeWidgetItem]): The list of items to search.
                text (str): The text to search for.
                field_index (int): The field index to search for.

            Yields:
                Generator[list[QTreeWidgetItem], None, None]: The list of nodes
                    that match the search criteria.

            """

            children: list[QTreeWidgetItem] = []

            for item in items:
                if (
                    item.text(field_index).__str__().lower().find(text.lower())
                    > -1
                ):
                    children.append(item)

                for child in [item.child(i) for i in range(item.childCount())]:
                    if (
                        child.text(field_index)
                        .__str__()
                        .lower()
                        .find(text.lower())
                        > -1
                    ):
                        children.append(child)

                    if child.childCount() > 0:
                        yield from recursive_find_in_children(
                            [child.child(i) for i in range(child.childCount())],
                            text,
                            field_index,
                        )

            yield children

        for item in [
            *recursive_find_in_children(
                items=self.__nodes_tree_widget.findItems(
                    None,
                    QtCore.Qt.MatchRecursive,
                    names_field_index,
                ),
                text=text,
                field_index=names_field_index,
            ),
            *recursive_find_in_children(
                items=self.__nodes_tree_widget.findItems(
                    None,
                    QtCore.Qt.MatchContains,
                    ids_field_index,
                ),
                text=text,
                field_index=ids_field_index,
            ),
        ]:
            name_items.extend(item)

        if len(name_items) > 0:
            self.__nodes_tree_widget.scrollToItem(
                name_items[0],
                QAbstractItemView.PositionAtCenter,
            )

        for item in set(name_items):
            self.__nodes_tree_widget.setItemSelected(item, True)

    def __build_main_menu(self, parent: Any) -> None:
        self.menubar = self.menuBar()
        self.file_menu = QMenu("File", parent)
        self.file_menu.addAction(
            "Import Tree", self.__load_tree_from_widget_menu
        ).setShortcut("Ctrl+O")

        self.menubar.addMenu(self.file_menu)

    def __load_tree_from_widget_menu(
        self,
    ) -> None:
        options = QFileDialog.Options()
        tree_filter = "Tree files (*.phylo.json);;All files (*.*)"

        qt_filename, _ = QFileDialog.getOpenFileName(
            self,
            "Open tree file",
            "",
            filter=tree_filter,
            selectedFilter=tree_filter,
            options=options,
        )

        if not qt_filename:
            return

        if (
            load_tree_either := self.__load_tree_from_phylojson_file(
                Path(qt_filename)
            )
        ).is_left:
            raise Exception(load_tree_either.value.msg)

        return

    def __load_tree_from_phylojson_file(
        self,
        tree_path: Path,
    ) -> Either[c_exc.MappedErrors, Literal[True]]:
        if not tree_path.exists():
            return c_exc.DadaTransferObjectError(
                f"Invalid path: {tree_path}",
                logger=LOGGER,
            )()

        if not tree_path.is_file():
            return c_exc.DadaTransferObjectError(
                f"Invalid path: {tree_path}",
                logger=LOGGER,
            )()

        if not tree_path.__str__().endswith(
            ("." + TreeSourceFormatEnum.PHYLO_JSON.value)
        ):
            return c_exc.DadaTransferObjectError(
                f"Invalid file type: {tree_path.suffix}",
                logger=LOGGER,
            )()

        with open(tree_path, "r") as f:
            self.__set_tree_widget_content(
                tree=ExtendedBioPythonTree.from_dict(content=load(f))
            )

        return right(True)

    def __set_tree_widget_content(
        self,
        tree: ExtendedBioPythonTree,
    ) -> None:
        root_node: ExtendedBioPythonClade
        self.__tree = tree
        self.__nodes_tree_widget.clear()

        def build_root_clade(
            clade: ExtendedBioPythonClade,
            parent: QTreeWidgetItem,
        ) -> QTreeWidgetItem:
            """Build the root node of the tree.

            Args:
                clade (ExtendedBioPythonClade): The root node.
                parent (QTreeWidgetItem): The parent node.

            Returns:
                QTreeWidgetItem: The root node.

            """

            parent.setExpanded(True)

            if clade.is_terminal():
                parent.setText(0, clade.name or "Unnamed clade")
            else:
                button = QPushButton()
                button.setText(clade.name or "Click to update")

                if clade.name is not None and clade.name != "":
                    button.setStyleSheet(
                        "QPushButton {background-color: #4CAF50;}"
                    )

                button.setFocusPolicy(QtCore.Qt.NoFocus)

                button.mouseDoubleClickEvent = (
                    lambda event: self.__annotate_node(
                        clade_id=clade._id,
                        current_value=clade.name,
                        event=event,
                    )
                )

                button.keyPressEvent = lambda event: self.__annotate_node(
                    clade_id=clade._id,
                    current_value=clade.name,
                    event=event,
                )

                parent.treeWidget().setItemWidget(parent, 0, button)

            for index, element in enumerate(
                [
                    (
                        "Clade Name or ID",
                        parent.setText,
                        ("" if clade.name is None else clade.name),
                    ),
                    (
                        "Clade Support",
                        parent.setText,
                        (
                            clade.confidence.__str__()
                            if clade.confidence
                            else ""
                        ),
                    ),
                    (
                        "If Clade is Outgroup",
                        parent.setText,
                        (
                            clade._is_outgroup.__str__()
                            if clade._is_outgroup
                            else ""
                        ),
                    ),
                    (
                        "If Clade is Terminal",
                        parent.setText,
                        (
                            clade.is_terminal().__str__()
                            if clade.is_terminal()
                            else ""
                        ),
                    ),
                    (
                        "The node unique identifier",
                        parent.setText,
                        (
                            clade._id.__str__()
                            if not clade.is_terminal()
                            else ""
                        ),
                    ),
                ]
            ):
                tip, action, value = element
                action(index, value)
                parent.setToolTip(index, tip)

            return parent

        def build_recursive_children(
            clade: ExtendedBioPythonClade,
            parent: QTreeWidgetItem,
        ) -> None:
            """Recursively build the tree from the root node.

            Args:
                clade (ExtendedBioPythonClade): The root node.
                parent (QTreeWidgetItem): The parent node.

            """

            child: ExtendedBioPythonClade
            for child in sorted(
                clade.clades,
                key=lambda x: x.clades.__len__(),
            ):
                item = QTreeWidgetItem(parent)
                child_item = build_root_clade(clade=child, parent=item)
                build_recursive_children(child, child_item)

            return

        root_item = QTreeWidgetItem(self.__nodes_tree_widget)
        root_node = self.__tree.root
        root_item = build_root_clade(clade=root_node, parent=root_item)

        build_recursive_children(root_node, root_item)

        return

    def __annotate_node(
        self,
        clade_id: UUID,
        current_value: str | None,
        event: QtGui.QMouseEvent,
    ) -> None:
        text, ok = QInputDialog.getText(
            self,
            "Annotate Node",
            "Enter node name:",
            QLineEdit.Normal,
            current_value or "",
        )

        if ok:
            self.__tree.set_clade_name(  # type: ignore
                clade_id=clade_id,
                name=(text if text != "" and text is not None else None),
            )

            if self.__tree is not None:
                with self.__tree_file_path.open("w") as f:
                    dump(self.__tree.to_dict(), f, indent=4, default=str)  # type: ignore

                with self.__tree_file_path.open("r") as f:
                    self.__set_tree_widget_content(
                        tree=ExtendedBioPythonTree.from_dict(content=load(f))
                    )

            item = self.__nodes_tree_widget.findItems(
                clade_id.__str__(), QtCore.Qt.MatchRecursive, 4
            )[0]

            self.__nodes_tree_widget.setItemSelected(item, True)
            self.__nodes_tree_widget.scrollToItem(
                item,
                QAbstractItemView.PositionAtCenter,
            )

        return


def main(phylo_json_tree: Path) -> None:
    """Main function.

    Args:
        phylo_json_tree (Path): The path to the phylo.json file.

    """

    app = QApplication(sys.argv)

    QtCore.QDir.addSearchPath(
        "assets",
        Path(__file__).absolute().parent.joinpath("assets").__str__(),
    )

    window = TreeEditor(tree_file_path=phylo_json_tree)
    window.show()
    app.setStyleSheet(STYLE_SHEET_DARK.read_text())
    sys.exit(app.exec_())
