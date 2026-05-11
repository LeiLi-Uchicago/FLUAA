#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from flu_pipeline.fasta import parse_isolate_from_seq_name
from flu_pipeline.nextclade import HA_COLUMNS_BY_SUBTYPE, gene_from_nextclade_dir_name, read_nextclade_tsv, resolve_column, seq_name_column


def main() -> None:
    args = parse_args()
    metadata = pd.read_csv(args.metadata, dtype=str, keep_default_na=False)
    id_column = "Isolate_Id" if "Isolate_Id" in metadata.columns else metadata.columns[0]

    for column in HA_COLUMNS_BY_SUBTYPE.get(args.subtype, []):
        metadata[column] = ""
    metadata["NA_clade"] = ""

    for directory in args.nextclade_dirs:
        gene = gene_from_nextclade_dir_name(Path(directory).name)
        tsv = Path(directory) / "nextclade.tsv"
        frame = read_nextclade_tsv(tsv)
        if frame.empty or gene not in {"HA", "NA"}:
            continue
        name_col = seq_name_column(frame)
        frame["_Isolate_Id"] = frame[name_col].map(parse_isolate_from_seq_name)
        frame = frame.drop_duplicates(subset=["_Isolate_Id"], keep="first")
        if gene == "NA":
            source_col = resolve_column(frame, "clade")
            if source_col:
                frame[source_col] = frame[source_col].fillna("")
                metadata = metadata.merge(
                    frame[["_Isolate_Id", source_col]].rename(columns={"_Isolate_Id": id_column, source_col: "_NA_clade"}),
                    on=id_column,
                    how="left",
                )
                metadata["NA_clade"] = metadata["_NA_clade"].fillna(metadata["NA_clade"])
                metadata = metadata.drop(columns=["_NA_clade"])
            continue

        for target in HA_COLUMNS_BY_SUBTYPE.get(args.subtype, []):
            logical = target.replace("HA_", "")
            source_col = resolve_column(frame, logical)
            if source_col is None:
                continue
            frame[source_col] = frame[source_col].fillna("")
            metadata = metadata.merge(
                frame[["_Isolate_Id", source_col]].rename(columns={"_Isolate_Id": id_column, source_col: f"_{target}"}),
                on=id_column,
                how="left",
            )
            metadata[target] = metadata[f"_{target}"].fillna(metadata[target])
            metadata = metadata.drop(columns=[f"_{target}"])

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(args.out, index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Annotate merged flu metadata with HA/NA Nextclade clades.")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--subtype", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--nextclade-dirs", nargs="+", required=True)
    return parser.parse_args()


if __name__ == "__main__":
    main()
