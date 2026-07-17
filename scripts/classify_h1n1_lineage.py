#!/usr/bin/env python3
"""Assign each H1N1 isolate to the seasonal or pdm09 lineage from its HA sequence.

The HA gene FASTA is aligned by Nextclade against both the pdm09 reference and the
custom seasonal reference (upstream, by the pipeline). This script reads the two
resulting Nextclade TSVs and assigns each isolate to whichever reference it aligns
to with the least divergence (substitutions + deletions + insertions). The lineage
then drives which reference set every one of that isolate's segments is mapped to.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from flu_pipeline.fasta import parse_isolate_from_seq_name
from flu_pipeline.nextclade import read_nextclade_tsv, seq_name_column


DIVERGENCE_COLUMNS = ["totalSubstitutions", "totalDeletions", "totalInsertions"]


def divergence_by_isolate(tsv_path: Path) -> dict[str, float]:
    """Map Isolate_Id -> divergence from the reference used for this Nextclade run.

    Lower means a closer match to that reference. A sequence that Nextclade could
    not align to this reference is scored as ``inf`` (maximally divergent), NOT 0.
    When a seasonal strain fails to align to the pdm09 reference, Nextclade still
    emits a row but leaves the mutation columns (and qc status) empty; summing
    those empties would wrongly read as a perfect match. An empty/failed alignment
    is itself evidence the isolate belongs to the OTHER reference's lineage, so it
    must never win the comparison.
    """
    frame = read_nextclade_tsv(tsv_path)
    if frame.empty:
        return {}
    name_col = seq_name_column(frame)
    status_col = "qc.overallStatus" if "qc.overallStatus" in frame.columns else None

    divergence: dict[str, float] = {}
    for row in frame.to_dict(orient="records"):
        isolate_id = parse_isolate_from_seq_name(str(row.get(name_col, "")))
        if not isolate_id:
            continue

        # Sum the divergence components, tracking whether any real value was present.
        total = 0.0
        have_data = False
        for column in DIVERGENCE_COLUMNS:
            value = str(row.get(column, "") or "").strip()
            if value != "":
                total += float(value)
                have_data = True

        status = str(row.get(status_col, "")).lower() if status_col else ""
        if not have_data or status == "bad":
            # No usable alignment against this reference -> maximally divergent.
            score = math.inf
        else:
            score = total

        # Keep the best (lowest-divergence) record if an isolate appears twice.
        if isolate_id not in divergence or score < divergence[isolate_id]:
            divergence[isolate_id] = score
    return divergence


def classify(pdm: dict[str, float], seasonal: dict[str, float]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for isolate_id in sorted(set(pdm) | set(seasonal)):
        pdm_div = pdm.get(isolate_id, math.inf)
        sea_div = seasonal.get(isolate_id, math.inf)
        # Ties and the "aligned to neither" case default to pdm09, matching the
        # pipeline's historical default for unclassifiable isolates.
        lineage = "seasonal" if sea_div < pdm_div else "pdm09"

        both_aligned = pdm_div != math.inf and sea_div != math.inf
        both_failed = pdm_div == math.inf and sea_div == math.inf
        if both_failed:
            # Aligned to neither reference: cannot classify (defaults to pdm09).
            confidence = "LOW"
            gap: float = math.inf
        elif not both_aligned:
            # Aligned cleanly to one reference and not at all to the other: the
            # references are divergent enough that this is a strong lineage call.
            confidence = "HIGH"
            gap = math.inf
        else:
            gap = abs(pdm_div - sea_div)
            if gap < 10:
                confidence = "AMBIGUOUS"
            elif gap >= 50:
                confidence = "HIGH"
            else:
                confidence = "MODERATE"

        rows.append(
            {
                "Isolate_Id": isolate_id,
                "H1_lineage": lineage,
                "confidence": confidence,
                "pdm_divergence": _fmt_div(pdm_div),
                "seasonal_divergence": _fmt_div(sea_div),
                "gap_size": _fmt(gap) if both_aligned else "",
            }
        )
    return rows


def _fmt(value: float) -> str:
    return "" if value == math.inf else f"{value:g}"


def _fmt_div(value: float) -> str:
    """Format a divergence score; a failed/empty alignment is shown explicitly
    rather than as a blank or a misleading number."""
    return "no_alignment" if value == math.inf else f"{value:g}"


def main() -> None:
    args = parse_args()
    pdm = divergence_by_isolate(Path(args.pdm_tsv))
    seasonal = divergence_by_isolate(Path(args.seasonal_tsv))
    rows = classify(pdm, seasonal)

    # Create output directory for classification results
    out_path = Path(args.out)
    out_dir = out_path.parent / "h1n1_classification"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Write main lineage CSV to the specified location
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Isolate_Id",
                "H1_lineage",
                "confidence",
                "pdm_divergence",
                "seasonal_divergence",
                "gap_size",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    # Write detailed divergence report to the classification folder
    detailed_path = out_dir / "divergence_details.csv"
    with detailed_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Isolate_Id",
                "H1_lineage",
                "confidence",
                "pdm_divergence",
                "seasonal_divergence",
                "gap_size",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    # Write summary statistics
    summary_path = out_dir / "classification_summary.txt"
    counts = {"pdm09": 0, "seasonal": 0, "HIGH": 0, "MODERATE": 0, "AMBIGUOUS": 0, "LOW": 0}
    for row in rows:
        counts[row["H1_lineage"]] += 1
        counts[row["confidence"]] += 1

    with summary_path.open("w") as handle:
        handle.write("=== H1N1 Lineage Classification Summary ===\n\n")
        handle.write(f"Total sequences: {len(rows)}\n")
        handle.write(f"  pdm09: {counts['pdm09']}\n")
        handle.write(f"  seasonal: {counts['seasonal']}\n\n")
        handle.write("Confidence distribution:\n")
        handle.write(f"  HIGH: {counts['HIGH']} sequences\n")
        handle.write(f"  MODERATE: {counts['MODERATE']} sequences\n")
        handle.write(f"  AMBIGUOUS: {counts['AMBIGUOUS']} sequences\n")
        handle.write(f"  LOW: {counts['LOW']} sequences\n\n")
        handle.write("References:\n")
        handle.write(f"  pdm09: {args.pdm_tsv}\n")
        handle.write(f"  seasonal: {args.seasonal_tsv}\n\n")
        handle.write("Output files:\n")
        handle.write(f"  {out_path.name}: Main lineage assignments (root level)\n")
        handle.write(f"  h1n1_classification/divergence_details.csv: Detailed divergence scores\n")
        handle.write(f"  h1n1_classification/classification_summary.txt: This summary\n")

    print(
        f"[classify_h1n1_lineage] pdm09={counts['pdm09']} seasonal={counts['seasonal']} "
        f"(HIGH={counts['HIGH']} MODERATE={counts['MODERATE']} AMBIGUOUS={counts['AMBIGUOUS']} LOW={counts['LOW']})"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Classify H1N1 isolates as seasonal or pdm09 from HA alignment divergence.")
    parser.add_argument("--pdm-tsv", required=True, help="Nextclade TSV of HA aligned to the pdm09 reference")
    parser.add_argument("--seasonal-tsv", required=True, help="Nextclade TSV of HA aligned to the seasonal reference")
    parser.add_argument("--out", required=True, help="Output lineage assignment CSV")
    return parser.parse_args()


if __name__ == "__main__":
    main()
