#!/usr/bin/env python3
"""Reorganize count and metadata outputs by H1N1 lineage into separate directories."""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    lineage_csv = Path(args.lineage_csv)
    metadata = Path(args.metadata)
    count_dir = Path(args.count_dir)

    # Read lineage assignments
    lineage_df = pd.read_csv(lineage_csv, dtype=str, keep_default_na=False)
    lineage_by_id = dict(zip(lineage_df["Isolate_Id"], lineage_df["H1_lineage"]))

    # Split metadata by lineage
    meta_df = pd.read_csv(metadata, dtype=str, keep_default_na=False)
    id_column = "Isolate_Id" if "Isolate_Id" in meta_df.columns else meta_df.columns[0]

    pdm_meta = meta_df[meta_df[id_column].apply(lambda x: lineage_by_id.get(str(x).strip(), "pdm09") == "pdm09")]
    seasonal_meta = meta_df[meta_df[id_column].apply(lambda x: lineage_by_id.get(str(x).strip(), "pdm09") == "seasonal")]

    # Create output directories
    pdm_dir = outdir / "pdm09"
    seasonal_dir = outdir / "seasonalH1N1"
    pdm_dir.mkdir(parents=True, exist_ok=True)
    seasonal_dir.mkdir(parents=True, exist_ok=True)

    # Write split metadata
    pdm_meta.to_csv(pdm_dir / "metadata_merged_annotated.csv", index=False)
    seasonal_meta.to_csv(seasonal_dir / "metadata_merged_annotated.csv", index=False)

    # Copy count directories, removing lineage suffix from protein names
    copy_lineage_counts(count_dir, pdm_dir / "count", "pdm09")
    copy_lineage_counts(count_dir, seasonal_dir / "count", "seasonal")

    print(f"[organize_by_lineage] Created pdm09/ and seasonalH1N1/ with split metadata and counts")


def copy_lineage_counts(src_count_dir: Path, dst_count_dir: Path, lineage: str) -> None:
    """Copy count files for a specific lineage, renaming protein directories by removing suffix."""
    src_count_dir = Path(src_count_dir)
    dst_count_dir = Path(dst_count_dir)

    # Find all protein directories with this lineage suffix
    for protein_dir in src_count_dir.glob(f"*_{lineage}"):
        base_protein = protein_dir.name.replace(f"_{lineage}", "")
        dst_protein_dir = dst_count_dir / base_protein
        dst_protein_dir.mkdir(parents=True, exist_ok=True)

        # Copy all CSV files from this protein's count directory
        for csv_file in protein_dir.glob("*.csv"):
            shutil.copy2(csv_file, dst_protein_dir / csv_file.name)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Reorganize outputs by H1N1 lineage.")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--lineage-csv", required=True)
    parser.add_argument("--count-dir", required=True)
    parser.add_argument("--outdir", default=".")
    return parser.parse_args()


if __name__ == "__main__":
    main()
