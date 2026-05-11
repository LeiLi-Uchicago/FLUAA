#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

from flu_pipeline.fasta import ambiguous_fraction, iter_fasta, ungapped_len, write_fasta_record
from flu_pipeline.metadata import read_metadata_files


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    metadata_dir = outdir / "metadata"
    fasta_dir = outdir / "fasta_by_gene"
    report_dir = outdir / "reports"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    fasta_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    metadata = read_metadata_files(Path(args.input_dir), args.id_column, args.date_column)
    metadata.to_csv(metadata_dir / "merged_metadata.csv", index=False)

    best: dict[tuple[str, str], dict[str, object]] = {}
    stats = Counter()
    conflict_rows: list[dict[str, object]] = []
    malformed_rows: list[dict[str, object]] = []
    processed = 0

    for fasta_path in sorted(path for path in Path(args.input_dir).glob("*-NT.fasta") if not path.name.startswith(".")):
        try:
            for record in iter_fasta(fasta_path):
                processed += 1
                if args.max_records_test and processed > args.max_records_test:
                    break
                stats["records_seen"] += 1
                key = (record.isolate_id, record.gene)
                candidate = {
                    "header": record.header,
                    "sequence": record.sequence,
                    "source_file": fasta_path.name,
                    "ungapped_len": ungapped_len(record.sequence),
                    "ambiguous_fraction": ambiguous_fraction(record.sequence),
                }
                previous = best.get(key)
                if previous is None:
                    best[key] = candidate
                    continue
                if previous["sequence"] == record.sequence:
                    stats["exact_duplicates"] += 1
                    continue
                stats["nonidentical_duplicates"] += 1
                chosen = choose_best(previous, candidate)
                best[key] = chosen
                conflict_rows.append(
                    {
                        "isolate_id": record.isolate_id,
                        "gene": record.gene,
                        "previous_header": previous["header"],
                        "new_header": candidate["header"],
                        "kept_header": chosen["header"],
                        "previous_ungapped_len": previous["ungapped_len"],
                        "new_ungapped_len": candidate["ungapped_len"],
                        "previous_ambiguous_fraction": previous["ambiguous_fraction"],
                        "new_ambiguous_fraction": candidate["ambiguous_fraction"],
                    }
                )
        except ValueError as exc:
            malformed_rows.append({"source_file": fasta_path.name, "error": str(exc)})
        if args.max_records_test and processed > args.max_records_test:
            break

    handles = {}
    try:
        for (_isolate_id, gene), row in sorted(best.items(), key=lambda item: (item[0][1], item[0][0])):
            handle = handles.get(gene)
            if handle is None:
                handle = (fasta_dir / f"{gene}.fasta").open("w")
                handles[gene] = handle
            write_fasta_record(handle, str(row["header"]), str(row["sequence"]))
    finally:
        for handle in handles.values():
            handle.close()

    write_dicts(report_dir / "fasta_duplicate_conflicts.csv", conflict_rows)
    write_dicts(report_dir / "malformed_fasta_headers.csv", malformed_rows)
    with (report_dir / "prepare_summary.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        for key, value in sorted(stats.items()):
            writer.writerow([key, value])
        writer.writerow(["unique_isolate_gene_records", len(best)])
        writer.writerow(["metadata_rows", len(metadata)])


def choose_best(previous: dict[str, object], candidate: dict[str, object]) -> dict[str, object]:
    previous_score = (int(previous["ungapped_len"]), -float(previous["ambiguous_fraction"]))
    candidate_score = (int(candidate["ungapped_len"]), -float(candidate["ambiguous_fraction"]))
    return candidate if candidate_score > previous_score else previous


def write_dicts(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge metadata and split/deduplicate flu FASTA records by gene.")
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--subtype", required=True)
    parser.add_argument("--id-column", default="Isolate_Id")
    parser.add_argument("--date-column", default="Collection_Date")
    parser.add_argument("--dedupe-mode", choices=["best_sequence"], default="best_sequence")
    parser.add_argument("--max-records-test", type=int, default=0)
    return parser.parse_args()


if __name__ == "__main__":
    main()
