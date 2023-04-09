#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Manual validation of the bayes algorithm when applied to the phylogenetic
tree.
"""


from json import load
from pathlib import Path
from uuid import UUID

import pandas as pd

from classeq.core.domain.dtos.priors import PriorGroup  # noqa: F401
from classeq.core.domain.dtos.clade import NodeType  # noqa: F401


def __load_kmers_map_and_linear_tree_as_dfs(
    reference_set_path: Path,
    train_source_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, list[tuple[UUID, pd.DataFrame]]]:
    """Load the kmer map and linear tree as pandas DataFrames.

    Args:
        reference_set_data (Path): Path to the reference set data.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: The kmer map and linear tree as
            pandas DataFrames.

    """
    kmer_map_df: pd.DataFrame
    linear_tree_df: pd.DataFrame
    internal_nodes_priors_dfs: list[tuple[UUID, pd.DataFrame]] = []

    with reference_set_path.open() as f:
        binding_data = load(f)

        labels_map = binding_data.get("labels_map")
        labels_map_df = pd.DataFrame([(k, v) for k, v in labels_map.items()])
        labels_map_df.columns = ["label", "index"]
        labels_map_df.set_index("index", inplace=True)

        binding_indices: dict[str, dict[int, int]] = dict()
        for index in (
            binding_data.get("msa").get("kmers_indices").get("indices")
        ):
            index_df = pd.DataFrame([1 for _ in index.get("records")])
            index_df.columns = [index.get("kmer")]
            index_df.index = index.get("records")
            binding_indices.update(index_df.to_dict())

        kmer_map_df = labels_map_df.join(
            pd.DataFrame(binding_indices).fillna(0).astype(int)
        ).sort_index()

        linear_tree_df = pd.DataFrame(binding_data.get("linear_tree"))
        linear_tree_df.set_index("id", inplace=True)
        linear_tree_df["type"] = linear_tree_df["type"].apply(lambda x: eval(x))

    with train_source_path.open() as f:
        binding_data = load(f)
        ingroups = binding_data.get("ingroups")

        for item in ingroups:
            parent: UUID
            if (parent := item.get("parent")) is None:
                raise ValueError("No parent found.")
            parent = UUID(parent)

            if (priors := item.get("priors")) is None:
                raise ValueError("No priors found.")

            try:
                ingroup_labels = next(
                    prior
                    for prior in priors
                    if (
                        "group" in prior
                        and eval(prior.get("group")) == PriorGroup.INGROUP
                    )
                )
            except StopIteration:
                raise ValueError("No ingroup priors found.")

            if (labels := ingroup_labels.get("labels")) is None:
                raise ValueError("No ingroup labels found.")

            ingroup_labels_kmers_source: dict[str, list[int]] = dict(
                {
                    "group": [PriorGroup.INGROUP for _ in labels],
                    "labels": labels,
                    **{
                        kmer: [1 for _ in labels]
                        for kmer, _ in ingroup_labels.get("priors").items()
                    },
                }
            )

            try:
                sister_labels = next(
                    prior
                    for prior in priors
                    if (
                        "group" in prior
                        and eval(prior.get("group")) == PriorGroup.SISTER
                    )
                )
            except StopIteration:
                raise ValueError("No sister-group priors found.")

            if (labels := sister_labels.get("labels")) is None:
                raise ValueError("No ingroup labels found.")

            sister_labels_kmers_source: dict[str, list[int]] = dict(
                {
                    "group": [PriorGroup.SISTER for _ in labels],
                    "labels": labels,
                    **{
                        kmer: [1 for _ in labels]
                        for kmer, _ in sister_labels.get("priors").items()
                    },
                }
            )

            item_priors = pd.concat(
                [
                    pd.DataFrame(ingroup_labels_kmers_source),
                    pd.DataFrame(sister_labels_kmers_source),
                ],
                axis=0,
            )

            internal_nodes_priors_dfs.append((parent, item_priors))

    return (kmer_map_df, linear_tree_df, internal_nodes_priors_dfs)


if __name__ == "__main__":
    (
        kmers_map,
        linear_tree,
        internal_nodes_priors_dfs,
    ) = __load_kmers_map_and_linear_tree_as_dfs(
        reference_set_path=Path(
            "Colletotrichum_acutatum_gapdh.sanitized.reference-set.json"
        ),
        train_source_path=Path(
            "Colletotrichum_acutatum_gapdh.sanitized.train-source.json"
        ),
    )

    kmers_map.to_csv("kmers_map.tsv", sep="\t")
    linear_tree.to_csv("linear_tree.tsv", sep="\t")

    for parent, priors_df in internal_nodes_priors_dfs:
        priors_df.to_csv(
            Path("raw-dfs").joinpath(parent.hex).with_suffix(".tsv"),
            sep="\t",
            index=False,
        )

        break
