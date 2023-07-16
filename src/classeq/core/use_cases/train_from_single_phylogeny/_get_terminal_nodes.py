from typing import Iterator

from classeq.core.domain.dtos.clade import CladeWrapper


def get_terminal_nodes(
    target_nodes: list[CladeWrapper],
    reference_nodes: list[CladeWrapper],
) -> list[CladeWrapper]:
    """Recursively get nodes of type NodeType.TERMINAL.

    Args:
        target_nodes (list[CladeWrapper]): A list of clades to find terminals
            from.
        reference_nodes (list[CladeWrapper]): A list of the remaining nodes to
            find complementary terminal nodes.

    Returns:
        list[CladeWrapper]: A list of clades of the type NodeType.TERMINAL.

    """

    def get_children(
        target_nodes: list[CladeWrapper],
    ) -> Iterator[CladeWrapper]:
        """Performs a recursive search for terminal nodes.

        Args:
            target_nodes (list[CladeWrapper]): A list of the clades to search
                terminals from.

        Yields:
            Iterator[CladeWrapper]: A single clade of the type
                NodeType.TERMINAL.
        """

        node: CladeWrapper
        for node in target_nodes:
            if node.is_terminal() or node.is_outgroup():
                yield node

            yield from get_terminal_nodes(
                target_nodes=[
                    i for i in reference_nodes if i.parent == node.id
                ],
                reference_nodes=[i for i in reference_nodes if i.id != node.id],
            )

    return [i for i in get_children(target_nodes)]
