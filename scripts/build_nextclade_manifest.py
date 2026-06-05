#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import pandas as pd

from flu_pipeline.fasta import iter_fasta, write_fasta_record
from flu_pipeline.metadata import parse_na_subtype


GENES = ["HA", "NA", "MP", "NP", "NS", "PA", "PB1", "PB2"]


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    split_dir = outdir / "split_fastas"
    split_dir.mkdir(parents=True, exist_ok=True)

    datasets = json.loads(args.datasets_json)
    seasonal_datasets = json.loads(args.h1n1_seasonal_datasets_json or "{}")
    h5_na_datasets = json.loads(args.h5_na_nextclade_datasets_json or "{}")
    metadata = pd.read_csv(args.metadata, dtype=str, keep_default_na=False)
    group_by_id = lineage_group_by_isolate(metadata, args.id_column, args.subtype)
    h5_na_by_id = h5_na_subtype_by_isolate(metadata, args.id_column)

    manifest_rows: list[dict[str, str]] = []
    for fasta in sorted(Path(path) for path in args.gene_fastas):
        gene = fasta.stem.upper()
        if gene not in GENES:
            continue
        if args.subtype == "H1N1" and args.h1n1_split_lineage:
            rows = split_h1n1_gene(fasta, split_dir, gene, group_by_id, datasets, seasonal_datasets)
        elif args.subtype == "H5NX" and gene == "NA":
            rows = split_h5_na_gene(fasta, split_dir, h5_na_by_id, h5_na_datasets)
        else:
            rows = write_unsplit_gene(fasta, split_dir, gene, datasets)
        manifest_rows.extend(rows)

    manifest_path = outdir / "nextclade_manifest.csv"
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["group", "gene", "fasta", "dataset"])
        writer.writeheader()
        writer.writerows(sorted(manifest_rows, key=lambda row: (row["group"], row["gene"])))


def lineage_group_by_isolate(metadata: pd.DataFrame, id_column: str, subtype: str) -> dict[str, str]:
    if subtype != "H1N1":
        return {}
    if id_column not in metadata.columns:
        id_column = "Isolate_Id" if "Isolate_Id" in metadata.columns else metadata.columns[0]
    groups: dict[str, str] = {}
    for row in metadata.to_dict(orient="records"):
        isolate_id = str(row.get(id_column, "")).strip()
        if not isolate_id:
            continue
        groups[isolate_id] = classify_h1n1_lineage(row)
    return groups


def classify_h1n1_lineage(row: dict[str, object]) -> str:
    lineage = str(row.get("Lineage", "") or "").strip().lower()
    year = parse_year(row.get("Year") or row.get("Collection_Date") or row.get("Collection date"))
    if lineage == "seasonal":
        return "seasonal"
    if lineage in {"pdm09", "pdm", "pandemic", "h1n1pdm09"}:
        return "pdm09"
    if year is not None and year < 2010:
        return "seasonal"
    return "pdm09"


def h5_na_subtype_by_isolate(metadata: pd.DataFrame, id_column: str) -> dict[str, str]:
    if id_column not in metadata.columns:
        id_column = "Isolate_Id" if "Isolate_Id" in metadata.columns else metadata.columns[0]
    groups: dict[str, str] = {}
    for row in metadata.to_dict(orient="records"):
        isolate_id = str(row.get(id_column, "")).strip()
        if not isolate_id:
            continue
        groups[isolate_id] = str(row.get("NA_subtype") or parse_na_subtype(row.get("Subtype", "")))
    return groups


def parse_year(value: object) -> int | None:
    text = str(value or "").strip()
    if len(text) >= 4 and text[:4].isdigit():
        return int(text[:4])
    return None


def split_h1n1_gene(
    fasta: Path,
    split_dir: Path,
    gene: str,
    group_by_id: dict[str, str],
    pdm09_datasets: dict[str, str],
    seasonal_datasets: dict[str, str],
) -> list[dict[str, str]]:
    grouped: dict[str, list[tuple[str, str]]] = {"pdm09": [], "seasonal": []}
    for record in iter_fasta(fasta):
        group = group_by_id.get(record.isolate_id, "pdm09")
        grouped.setdefault(group, []).append((record.header, record.sequence))

    rows: list[dict[str, str]] = []
    for group, records in grouped.items():
        if not records:
            continue
        dataset = (seasonal_datasets if group == "seasonal" else pdm09_datasets).get(gene)
        if not dataset:
            continue
        path = split_dir / group / f"{gene}.fasta"
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w") as handle:
            for header, sequence in records:
                write_fasta_record(handle, header, sequence)
        rows.append({"group": group, "gene": gene, "fasta": str(path.resolve()), "dataset": str(dataset)})
    return rows


def split_h5_na_gene(
    fasta: Path,
    split_dir: Path,
    subtype_by_id: dict[str, str],
    h5_na_datasets: dict[str, str],
) -> list[dict[str, str]]:
    grouped: dict[str, list[tuple[str, str]]] = {}
    for record in iter_fasta(fasta):
        group = subtype_by_id.get(record.isolate_id, "unknown")
        grouped.setdefault(group, []).append((record.header, record.sequence))

    rows: list[dict[str, str]] = []
    for group, records in sorted(grouped.items()):
        dataset = h5_na_datasets.get(group)
        if not dataset:
            continue
        path = split_dir / group / "NA.fasta"
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w") as handle:
            for header, sequence in records:
                write_fasta_record(handle, header, sequence)
        rows.append({"group": group, "gene": "NA", "fasta": str(path.resolve()), "dataset": str(dataset)})
    return rows


def write_unsplit_gene(
    fasta: Path,
    split_dir: Path,
    gene: str,
    datasets: dict[str, str],
) -> list[dict[str, str]]:
    dataset = datasets.get(gene)
    if not dataset:
        return []
    path = split_dir / "all" / f"{gene}.fasta"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(fasta.read_text())
    return [{"group": "all", "gene": gene, "fasta": str(path.resolve()), "dataset": str(dataset)}]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Nextclade input manifest, optionally splitting H1N1 by lineage.")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--gene-fastas", nargs="+", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--subtype", required=True)
    parser.add_argument("--id-column", default="Isolate_Id")
    parser.add_argument("--datasets-json", required=True)
    parser.add_argument("--h1n1-seasonal-datasets-json", default="{}")
    parser.add_argument("--h5-na-nextclade-datasets-json", default="{}")
    parser.add_argument("--h1n1-split-lineage", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    main()
