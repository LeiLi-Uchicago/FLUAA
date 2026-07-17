#!/usr/bin/env python3
"""Split metadata and FASTA files by H1N1 lineage into pdm09 and seasonal groups."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from flu_pipeline.fasta import iter_fasta, write_fasta_record


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)

    # Read lineage assignments
    lineage_df = pd.read_csv(args.lineage_csv, dtype=str, keep_default_na=False)
    lineage_by_id = dict(zip(lineage_df["Isolate_Id"], lineage_df["H1_lineage"]))

    # Read metadata and split
    metadata = pd.read_csv(args.metadata, dtype=str, keep_default_na=False)
    id_column = "Isolate_Id" if "Isolate_Id" in metadata.columns else metadata.columns[0]

    pdm_meta = []
    seasonal_meta = []
    pdm_ids = set()
    seasonal_ids = set()

    for _, row in metadata.iterrows():
        isolate_id = str(row[id_column]).strip()
        lineage = lineage_by_id.get(isolate_id, "pdm09")  # Default to pdm09
        if lineage == "seasonal":
            seasonal_meta.append(row)
            seasonal_ids.add(isolate_id)
        else:
            pdm_meta.append(row)
            pdm_ids.add(isolate_id)

    # Write split metadata
    outdir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(pdm_meta).to_csv(outdir / "pdm09_metadata.csv", index=False)
    pd.DataFrame(seasonal_meta).to_csv(outdir / "seasonal_metadata.csv", index=False)

    print(f"[split_by_lineage] pdm09={len(pdm_ids)}, seasonal={len(seasonal_ids)}")

    # Split FASTA files
    pdm_fastas = outdir / "pdm09_fastas"
    seasonal_fastas = outdir / "seasonal_fastas"
    pdm_fastas.mkdir(exist_ok=True)
    seasonal_fastas.mkdir(exist_ok=True)

    for fasta_path in args.gene_fastas:
        fasta_path = Path(fasta_path)
        gene = fasta_path.stem

        pdm_out = pdm_fastas / f"{gene}.fasta"
        seasonal_out = seasonal_fastas / f"{gene}.fasta"

        pdm_handle = pdm_out.open("w")
        seasonal_handle = seasonal_out.open("w")

        try:
            for record in iter_fasta(fasta_path):
                isolate_id = record.isolate_id
                if isolate_id in seasonal_ids:
                    write_fasta_record(seasonal_handle, record.header, record.sequence)
                else:
                    write_fasta_record(pdm_handle, record.header, record.sequence)
        finally:
            pdm_handle.close()
            seasonal_handle.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Split metadata and FASTA by H1N1 lineage.")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--lineage-csv", required=True)
    parser.add_argument("--gene-fastas", nargs="+", required=True)
    parser.add_argument("--outdir", default=".")
    return parser.parse_args()


if __name__ == "__main__":
    main()
