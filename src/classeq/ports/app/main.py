import re
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
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from classeq.core.domain.dtos.biopython_wrappers import (
    ExtendedBioPythonClade,
    ExtendedBioPythonTree,
    ExtraTaxonomicRanks,
    MajorTaxonomicRanks,
    MinorTaxonomicRanks,
    try_to_reach_rank_enum,
)
from classeq.core.domain.dtos.tree_source_format import TreeSourceFormatEnum
from classeq.ports.app.settings import STYLE_SHEET_DARK
from classeq.settings import LOGGER


class TreeEditor(QMainWindow):
    # ? ------------------------------------------------------------------------
    # ? CLASS ATTRIBUTES
    # ? ------------------------------------------------------------------------

    __tree: ExtendedBioPythonTree | None = None
    __clade: ExtendedBioPythonClade | None = None
    __tree_file_path: Path
    __tree_field_names: int = 0
    __tree_field_status: int = 1
    __tree_field_taxid: int = 3
    __tree_field_related_rank: int = 4
    __tree_field_ids: int = 7

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

        central_item = QWidget(self)
        self.setCentralWidget(central_item)

        self.__add_widgets(parent=central_item)
        self.__build_main_menu(self)

        self.__tree_file_path = tree_file_path

        if (
            load_tree_either := self.__load_tree_from_phylojson_file(
                tree_file_path
            )
        ).is_left:
            raise Exception(load_tree_either.value.msg)

    # ? ------------------------------------------------------------------------
    # ? SUPER METHODS REPLACEMENT
    # ? ------------------------------------------------------------------------

    def closeEvent(self, event: Any) -> None:
        """Override the close event to save the tree before closing the
        application.

        """

        LOGGER.info("Closing application")
        event.accept()

    # ? ------------------------------------------------------------------------
    # ? PRIVATE STATIC METHODS
    # ? ------------------------------------------------------------------------

    def __add_widgets(self, parent: Any) -> None:
        """Add the widgets to the main window."""

        # ? --------------------------------------------------------------------
        # ? Main layout
        #
        # Here the phylogeny is presented.
        #
        # ? --------------------------------------------------------------------
        layout = QVBoxLayout(parent)
        splitter = QSplitter(parent)

        # ? --------------------------------------------------------------------
        # ? Right panel = Details + Annotations
        # ? --------------------------------------------------------------------

        left_splitter = QSplitter(splitter)
        left_splitter.setOrientation(QtCore.Qt.Vertical)

        left_widget = QWidget(left_splitter)
        left_layout = QVBoxLayout(left_widget)
        print(left_layout)
        # left_splitter.addWidget(left_widget)

        table_fields = [
            "Node Name",
            "Clade Support",
            "Taxid",
            "Related Rank",
            "Outgroup",
            "Terminal",
            "Identifier",
            "Children",
        ]

        table_view = self.__table_view = QTableWidget(
            len(table_fields), 1, left_widget
        )

        table_view.setHorizontalHeaderLabels(["Values"])
        table_view.setVerticalHeaderLabels(table_fields)
        table_view.horizontalHeader().setStretchLastSection(True)
        table_view.verticalHeader().setStretchLastSection(True)
        table_view.setFixedHeight(280)
        table_view.resizeColumnsToContents()

        line_edit = QLineEdit(left_widget)

        left_splitter.addWidget(self.__table_view)
        left_splitter.addWidget(line_edit)

        # ? --------------------------------------------------------------------
        # ? Left panel = Filter box + Tree box
        # ? --------------------------------------------------------------------

        right_widget = QWidget(splitter)
        right_layout = QVBoxLayout(right_widget)
        splitter.addWidget(right_widget)

        search_field = self.__search_field = QLineEdit(right_widget)
        search_field.setFocusPolicy(QtCore.Qt.StrongFocus)
        search_field.setObjectName("search_field")
        search_field.setStyleSheet("font-size: 16px; height: 50px;")
        search_field.setPlaceholderText("Search for a node")
        search_field.textChanged.connect(
            lambda text: self.__on_text_changed(text)
        )

        tree_widget = self.__tree_widget = QTreeWidget(right_widget)
        tree_widget.setHeaderItem(
            QTreeWidgetItem(
                [
                    "Name",
                    " ",
                    "Support",
                    "Taxid",
                    "Related Rank",
                    "Outgroup",
                    "Terminal",
                    "Identifier",
                ]
            )
        )

        tree_widget.header().swapSections(0, 1)
        tree_widget.setColumnCount(8)
        tree_widget.setColumnWidth(0, 700)
        tree_widget.setColumnWidth(1, 10)
        tree_widget.setColumnWidth(2, 80)
        tree_widget.setColumnWidth(3, 80)
        tree_widget.setColumnWidth(4, 110)
        tree_widget.setColumnWidth(5, 80)
        tree_widget.setColumnWidth(6, 80)
        tree_widget.setColumnWidth(7, 300)
        tree_widget.setSelectionMode(QAbstractItemView.ExtendedSelection)
        tree_widget.setSortingEnabled(True)
        tree_widget.itemClicked.connect(self.__on_item_clicked)

        right_layout.addWidget(self.__search_field)
        right_layout.addWidget(self.__tree_widget)

        # ? --------------------------------------------------------------------
        # ? Update splitter layout
        # ? --------------------------------------------------------------------

        splitter.setStretchFactor(0, 2)
        splitter.setStretchFactor(1, 5)
        layout.addWidget(splitter)

        return

    @QtCore.Slot(QTreeWidgetItem, int)
    def __on_item_clicked(self, item: QTreeWidgetItem, _: int) -> None:
        if self.__tree is None:
            return

        self.__clade = self.__tree.find_clade_by_id(
            id=UUID(item.text(self.__tree_field_ids))
        )

        self.__set_table_data()

    def __on_text_changed(self, text: str) -> None:
        self.__tree_widget.clearSelection()

        if text == "" or text is None or text.__len__() < 3:
            return

        name_items: list[QTreeWidgetItem] = []

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
                items=self.__tree_widget.findItems(
                    None,
                    QtCore.Qt.MatchRecursive,
                    self.__tree_field_names,
                ),
                text=text,
                field_index=self.__tree_field_names,
            ),
            *recursive_find_in_children(
                items=self.__tree_widget.findItems(
                    None,
                    QtCore.Qt.MatchContains,
                    self.__tree_field_ids,
                ),
                text=text,
                field_index=self.__tree_field_ids,
            ),
        ]:
            name_items.extend(item)

        if len(name_items) > 0:
            self.__tree_widget.scrollToItem(
                name_items[0],
                QAbstractItemView.PositionAtCenter,
            )

        for item in set(name_items):
            self.__tree_widget.setItemSelected(item, True)

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

    @staticmethod
    def __check_status(
        clade: ExtendedBioPythonClade,
    ) -> QtGui.QIcon | None:
        """Check if the clade is valid.

        Args:
            clade (ExtendedBioPythonClade): The clade to check.

        Returns:
            QtGui.QIcon | None: The icon to show.

        """

        if clade.is_terminal():
            return None

        if clade.name is None:
            return QtGui.QIcon("assets:/icons/warning-16x16.png")
        return QtGui.QIcon("assets:/icons/check-simple-16x16.png")

    def __set_tree_widget_content(
        self,
        tree: ExtendedBioPythonTree,
    ) -> None:
        root_node: ExtendedBioPythonClade
        self.__tree = tree
        self.__tree_widget.clear()

        LOGGER.info("Building tree widget")

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
            faded_color = QtGui.QColor("#808080")

            if clade.is_terminal():
                parent.setTextColor(self.__tree_field_names, faded_color)

            else:
                # ? ------------------------------------------------------------
                # ? Annotate name
                # ? ------------------------------------------------------------

                name_button = QPushButton()
                name_button.setFocusPolicy(QtCore.Qt.NoFocus)
                name_button.setCursor(
                    QtGui.QCursor(QtCore.Qt.PointingHandCursor)
                )

                name_button.mouseDoubleClickEvent = (
                    lambda _: self.__annotate_clade_name(
                        clade=clade,
                        current_value=clade.name,
                    )
                )

                parent.treeWidget().setItemWidget(
                    parent, self.__tree_field_names, name_button
                )

                # ? ------------------------------------------------------------
                # ? Annotate taxid
                # ? ------------------------------------------------------------

                taxid_button = QPushButton()
                taxid_button.setFocusPolicy(QtCore.Qt.NoFocus)
                taxid_button.setCursor(
                    QtGui.QCursor(QtCore.Qt.PointingHandCursor)
                )

                taxid_button.mouseDoubleClickEvent = (
                    lambda _: self.__annotate_clade_taxid(
                        clade=clade,
                        current_value=clade._taxid.__str__(),
                    )
                )

                parent.treeWidget().setItemWidget(
                    parent, self.__tree_field_taxid, taxid_button
                )

                # ? ------------------------------------------------------------
                # ? Annotate name
                # ? ------------------------------------------------------------

                rank_button = QPushButton()
                rank_button.setFocusPolicy(QtCore.Qt.NoFocus)
                rank_button.setCursor(
                    QtGui.QCursor(QtCore.Qt.PointingHandCursor)
                )

                rank_button.mouseDoubleClickEvent = (
                    lambda _: self.__annotate_clade_rank(
                        clade=clade,
                        current_value=clade._related_rank,
                    )
                )

                parent.treeWidget().setItemWidget(
                    parent, self.__tree_field_related_rank, rank_button
                )

            for index, element in enumerate(
                [
                    (
                        "Clade Name or ID",
                        parent.setText,
                        (
                            clade.name.__str__()
                            if clade.name is not None
                            else None
                        ),
                    ),
                    (
                        "The clade general status",
                        (
                            parent.setText
                            if self.__check_status(clade) is None
                            else parent.setIcon
                        ),
                        self.__check_status(clade),
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
                        "Clade Related Taxid",
                        parent.setText,
                        (clade._taxid.__str__() if clade._taxid else ""),
                    ),
                    (
                        "Clade Taxonomic Rank",
                        parent.setText,
                        (
                            clade._related_rank.name
                            if clade._related_rank
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
                        clade._id.__str__(),
                    ),
                ]
            ):
                tip, action, value = element
                action(index, value)
                parent.setTextColor(index, faded_color)
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

                if child.clades.__len__() > 0:
                    build_recursive_children(child, child_item)

            return

        root_item = QTreeWidgetItem(self.__tree_widget)
        root_node = self.__tree.root

        LOGGER.info("\tBuilding tree root")

        root_item = build_root_clade(clade=root_node, parent=root_item)

        LOGGER.info("\tBuilding tree children")

        build_recursive_children(root_node, root_item)

        LOGGER.info("\tSetting up tree detail view")

        self.__set_table_data(root_node)

        LOGGER.info("\tBuilding done")

        return

    def __set_table_data(
        self,
        clade: ExtendedBioPythonClade | None = None,
    ) -> None:
        if clade is None:
            if self.__clade is None:
                return
            clade = self.__clade

        for index, data in enumerate(
            [
                (clade.name.__str__() if clade.name else ""),
                (clade.confidence.__str__() if clade.confidence else ""),
                (clade._taxid.__str__() if clade._taxid else ""),
                (clade._related_rank.name if clade._related_rank else ""),
                (clade._is_outgroup.__str__() or ""),
                (clade.is_terminal().__str__() if clade.is_terminal() else ""),
                (clade._id.__str__() or ""),
                (
                    clade.clades.__len__().__str__()
                    if not clade.is_terminal()
                    else ""
                ),
            ]
        ):
            self.__table_view.setItem(index, 0, QTableWidgetItem(data))
            item = self.__table_view.item(index, 0)
            item.setFlags(QtCore.Qt.ItemIsEnabled)

    def __annotate_clade_name(
        self,
        clade: ExtendedBioPythonClade,
        current_value: str | None,
    ) -> None:
        clade_id: UUID = clade._id

        text, ok = QInputDialog.getText(
            self,
            "Annotate Node",
            "Enter node name:",
            QLineEdit.Normal,
            current_value or "",
        )

        if ok:
            slug_text = re.sub(r"[^a-zA-Z0-9\s_]", "", text).replace(" ", "_")

            LOGGER.info("Updating node name")
            self.__tree.set_clade_name(  # type: ignore
                id=clade_id,
                name=(
                    slug_text
                    if slug_text != "" and slug_text is not None
                    else None
                ),
            )

            LOGGER.info("Persisting tree")
            if self.__tree is not None:
                LOGGER.info("\tSaving tree")
                with self.__tree_file_path.open("w") as f:
                    dump(self.__tree.to_dict(), f, indent=4, default=str)  # type: ignore

            LOGGER.info("Persisting Finished")

            self.__finish_annotation(
                clade=clade,
                column=self.__tree_field_names,
                text=slug_text,
            )

            LOGGER.info("Done")

        return

    def __annotate_clade_taxid(
        self,
        clade: ExtendedBioPythonClade,
        current_value: str | None,
    ) -> None:
        clade_id: UUID = clade._id

        text, ok = QInputDialog.getText(
            self,
            "Annotate Node",
            "Enter node name:",
            QLineEdit.Normal,
            current_value or "",
        )

        if ok:
            try:
                slug_text = int(text)
            except ValueError:
                slug_text = None

            LOGGER.info("Updating node name")
            self.__tree.set_clade_taxid(  # type: ignore
                id=clade_id,
                taxid=(slug_text if slug_text is not None else None),
            )

            LOGGER.info("Persisting tree")
            if self.__tree is not None:
                LOGGER.info("\tSaving tree")
                with self.__tree_file_path.open("w") as f:
                    dump(self.__tree.to_dict(), f, indent=4, default=str)  # type: ignore

            LOGGER.info("Persisting Finished")

            self.__finish_annotation(
                clade=clade,
                column=self.__tree_field_taxid,
                text=slug_text,
            )

            LOGGER.info("Done")

        return

    def __annotate_clade_rank(
        self,
        clade: ExtendedBioPythonClade,
        current_value: str | None,
    ) -> None:
        clade_id: UUID = clade._id

        text, ok = QInputDialog.getItem(
            self,
            "Annotate Node",
            "Enter node name:",
            [
                ExtraTaxonomicRanks.NO_RANK.value,
                *[rank.value for rank in MajorTaxonomicRanks],
                *[rank.value for rank in MinorTaxonomicRanks],
                *[rank.value for rank in ExtraTaxonomicRanks],
            ],
        )

        if ok:
            rank = try_to_reach_rank_enum(text)

            LOGGER.info("Updating node name")
            self.__tree.set_clade_rank(  # type: ignore
                id=clade_id,
                rank=rank,
            )

            LOGGER.info("Persisting tree")
            if self.__tree is not None:
                LOGGER.info("\tSaving tree")
                with self.__tree_file_path.open("w") as f:
                    dump(self.__tree.to_dict(), f, indent=4, default=str)  # type: ignore

            LOGGER.info("Persisting Finished")

            self.__finish_annotation(
                clade=clade,
                column=self.__tree_field_related_rank,
                text=rank.name,
            )

            LOGGER.info("Done")

        return

    def __finish_annotation(
        self,
        clade: ExtendedBioPythonClade,
        column: int,
        text: str | int | None,
    ) -> None:
        item = self.__tree_widget.findItems(
            clade._id.__str__(),
            QtCore.Qt.MatchRecursive,
            self.__tree_field_ids,
        )[0]

        item.setText(column, text.__str__())
        item.setIcon(self.__tree_field_status, self.__check_status(clade))

        self.__tree_widget.clearSelection()
        self.__tree_widget.setItemSelected(item, True)
        self.__tree_widget.scrollToItem(
            item,
            QAbstractItemView.PositionAtCenter,
        )

        self.__set_table_data(clade=clade)


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
