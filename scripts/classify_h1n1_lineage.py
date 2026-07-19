#!/usr/bin/env python3
"""Assign each H1N1 isolate to the seasonal or pdm09 lineage from its HA sequence.

The HA gene FASTA is aligned by Nextclade against both the pdm09 reference and the
custom seasonal reference (upstream, by the pipeline). This script reads the two
resulting Nextclade TSVs and assigns each isolate to whichever reference it aligns
to with the least divergence (substitutions + deletions + insertions). The lineage
then drives which reference set every one of that isolate's segments is mapped to.

Two guards keep the call honest:

1. Absolute/normalized divergence gate. A win only counts if the winning
   alignment is genuinely close to its reference. pdm09 and seasonal HA differ by
   ~22% at the nucleotide level (~380 substitutions over a full HA), so a "win"
   with a per-site divergence near that distance is not a real match -- it is a
   pdm09 sequence measured against the seasonal reference (or vice versa) after
   the other reference's alignment failed. Divergence is normalized by the
   aligned length so partial sequences are judged on the same scale. When the
   winning side is above the threshold the isolate is not confidently
   classifiable and is filtered out (written to lineage_filter_out.csv).

2. Pre-pandemic filter. Seasonal H1N1 stopped circulating in humans after the
   2009 pandemic, so any isolate that lands in pdm09 but was collected before
   ``--pdm-min-year`` is filtered out (written to lineage_filter_out.csv and
   excluded from both lineage outputs).
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import pandas as pd

from flu_pipeline.fasta import parse_isolate_from_seq_name
from flu_pipeline.nextclade import read_nextclade_tsv, seq_name_column


DIVERGENCE_COLUMNS = ["totalSubstitutions", "totalDeletions", "totalInsertions"]

# A finite divergence score paired with the aligned length (in reference
# nucleotides) it was measured over. ``span == 0`` means the aligned length is
# unknown, so the caller must fall back to an absolute threshold.
DivergenceScore = tuple[float, float]


def divergence_by_isolate(tsv_path: Path) -> dict[str, DivergenceScore]:
    """Map Isolate_Id -> (divergence, aligned_length) from a Nextclade run.

    Lower divergence means a closer match to that reference. A sequence that
    Nextclade could not align to this reference is scored as ``inf`` (maximally
    divergent), NOT 0. When a strain fails to align to the wrong-lineage
    reference, Nextclade still emits a row but leaves the mutation columns empty;
    summing those empties would wrongly read as a perfect match. An empty/failed
    alignment is itself evidence the isolate belongs to the OTHER reference's
    lineage, so it must never win the comparison.

    ``qc.overallStatus`` is deliberately NOT used to gate the score. That status
    reflects sequence QUALITY (e.g. too many mixed/ambiguous bases), not whether
    the alignment succeeded: a genuine pdm09 sequence full of ambiguity codes
    aligns to the pdm09 reference at low divergence yet is flagged "bad". Scoring
    such a sequence as ``inf`` would discard a correct, low-divergence alignment
    and push it to the wrong lineage. The actual divergence (with the
    normalized-divergence gate in ``classify``) is the honest signal instead.
    """
    frame = read_nextclade_tsv(tsv_path)
    if frame.empty:
        return {}
    name_col = seq_name_column(frame)

    scores: dict[str, DivergenceScore] = {}
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

        if not have_data:
            # No alignment against this reference -> maximally divergent.
            score = math.inf
            span = 0.0
        else:
            score = total
            span = aligned_span(row)

        # Keep the best (lowest-divergence) record if an isolate appears twice.
        if isolate_id not in scores or score < scores[isolate_id][0]:
            scores[isolate_id] = (score, span)
    return scores


def aligned_span(row: dict[str, object]) -> float:
    """Aligned length in reference nucleotides (alignmentEnd - alignmentStart).

    Used to normalize divergence so a short fragment and a full-length sequence
    are compared on the same per-site scale. Returns 0 when unavailable.
    """
    start = _to_float(row.get("alignmentStart"))
    end = _to_float(row.get("alignmentEnd"))
    if start is None or end is None:
        return 0.0
    span = end - start
    return span if span > 0 else 0.0


def _to_float(value: object) -> float | None:
    text = str(value if value is not None else "").strip()
    if text == "":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def classify(
    pdm: dict[str, DivergenceScore],
    seasonal: dict[str, DivergenceScore],
    year_by_id: dict[str, int],
    max_divergence_frac: float,
    max_divergence_abs: float,
    pdm_min_year: int,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for isolate_id in sorted(set(pdm) | set(seasonal)):
        pdm_div, pdm_span = pdm.get(isolate_id, (math.inf, 0.0))
        sea_div, sea_span = seasonal.get(isolate_id, (math.inf, 0.0))

        # Ties and the "aligned to neither" case default to pdm09, matching the
        # pipeline's historical default for unclassifiable isolates.
        preferred = "seasonal" if sea_div < pdm_div else "pdm09"
        win_div, win_span = (sea_div, sea_span) if preferred == "seasonal" else (pdm_div, pdm_span)

        both_aligned = pdm_div != math.inf and sea_div != math.inf

        if winning_alignment_is_real(win_div, win_span, max_divergence_frac, max_divergence_abs):
            lineage = preferred
            gap: float = math.inf
            if not both_aligned:
                # Aligned cleanly to one reference and not at all to the other,
                # and that one alignment passed the divergence gate: a strong call.
                confidence = "HIGH"
            else:
                gap = abs(pdm_div - sea_div)
                if gap < 10:
                    confidence = "AMBIGUOUS"
                elif gap >= 50:
                    confidence = "HIGH"
                else:
                    confidence = "MODERATE"
            filter_reason = ""
            if lineage == "pdm09":
                year = year_by_id.get(isolate_id)
                if year is not None and year < pdm_min_year:
                    # Seasonal H1N1 vanished from humans after 2009; a pre-pandemic
                    # isolate cannot genuinely be pdm09. Filter it out.
                    lineage = "filtered_out"
                    filter_reason = f"pdm09 but collected before {pdm_min_year} (year={year})"
        else:
            # The best alignment is too divergent to be a real match to either
            # reference (e.g. a pdm09 sequence that failed the pdm09 alignment and
            # only "aligned" to seasonal at the ~22% pdm-vs-seasonal distance), or
            # it aligned to neither. Cannot classify; filter it out.
            lineage = "filtered_out"
            confidence = "UNCLASSIFIED"
            gap = math.inf
            filter_reason = "unclassifiable: HA matches neither reference below the divergence threshold"

        rows.append(
            {
                "Isolate_Id": isolate_id,
                "H1_lineage": lineage,
                "confidence": confidence,
                "pdm_divergence": _fmt_div(pdm_div),
                "seasonal_divergence": _fmt_div(sea_div),
                "pdm_divergence_frac": _fmt_frac(pdm_div, pdm_span),
                "seasonal_divergence_frac": _fmt_frac(sea_div, sea_span),
                "gap_size": _fmt(gap) if both_aligned else "",
                "filter_reason": filter_reason,
            }
        )
    return rows


def winning_alignment_is_real(
    win_div: float, win_span: float, max_divergence_frac: float, max_divergence_abs: float
) -> bool:
    """Whether the winning alignment is close enough to be a genuine lineage match."""
    if win_div == math.inf:
        return False
    if win_span > 0:
        return (win_div / win_span) <= max_divergence_frac
    # Aligned length unknown: fall back to an absolute substitution threshold.
    return win_div <= max_divergence_abs


def year_by_isolate(metadata_path: str | None) -> dict[str, int]:
    if not metadata_path:
        return {}
    path = Path(metadata_path)
    if not path.exists() or path.stat().st_size == 0:
        return {}
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    if "Isolate_Id" not in frame.columns:
        return {}
    years: dict[str, int] = {}
    for row in frame.to_dict(orient="records"):
        isolate_id = str(row.get("Isolate_Id", "")).strip()
        if not isolate_id:
            continue
        year = parse_year(row.get("Year")) or parse_year(row.get("Collection_Date"))
        if year is not None:
            years[isolate_id] = year
    return years


def parse_year(value: object) -> int | None:
    text = str(value if value is not None else "").strip()
    if len(text) >= 4 and text[:4].isdigit():
        return int(text[:4])
    return None


def _fmt(value: float) -> str:
    return "" if value == math.inf else f"{value:g}"


def _fmt_div(value: float) -> str:
    """Format a divergence score; a failed/empty alignment is shown explicitly
    rather than as a blank or a misleading number."""
    return "no_alignment" if value == math.inf else f"{value:g}"


def _fmt_frac(value: float, span: float) -> str:
    """Format a per-site (normalized) divergence: divergence / aligned length.

    ``no_alignment`` when the alignment failed, blank when the aligned length is
    unknown (so it cannot be normalized)."""
    if value == math.inf:
        return "no_alignment"
    if span <= 0:
        return ""
    return f"{value / span:.4f}"


LINEAGE_FIELDS = [
    "Isolate_Id",
    "H1_lineage",
    "confidence",
    "pdm_divergence",
    "seasonal_divergence",
    "pdm_divergence_frac",
    "seasonal_divergence_frac",
    "gap_size",
    "filter_reason",
]


def main() -> None:
    args = parse_args()
    pdm = divergence_by_isolate(Path(args.pdm_tsv))
    seasonal = divergence_by_isolate(Path(args.seasonal_tsv))
    year_by_id = year_by_isolate(args.metadata)
    rows = classify(
        pdm,
        seasonal,
        year_by_id,
        max_divergence_frac=args.max_divergence_frac,
        max_divergence_abs=args.max_divergence_abs,
        pdm_min_year=args.pdm_min_year,
    )

    out_path = Path(args.out)
    out_dir = out_path.parent / "h1n1_classification"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Main lineage CSV (root level) and a detailed copy in the classification folder.
    for target in (out_path, out_dir / "divergence_details.csv"):
        with target.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=LINEAGE_FIELDS)
            writer.writeheader()
            writer.writerows(rows)

    # Filtered-out isolates: full metadata (when available) plus why they were dropped.
    filtered_rows = [row for row in rows if row["H1_lineage"] == "filtered_out"]
    write_filter_out(Path(args.filter_out), filtered_rows, args.metadata)

    # Summary statistics.
    counts = {"pdm09": 0, "seasonal": 0, "filtered_out": 0}
    conf_counts = {"HIGH": 0, "MODERATE": 0, "AMBIGUOUS": 0, "LOW": 0, "UNCLASSIFIED": 0}
    filtered_unclassifiable = 0
    filtered_pre_pandemic = 0
    for row in rows:
        counts[row["H1_lineage"]] = counts.get(row["H1_lineage"], 0) + 1
        conf_counts[row["confidence"]] = conf_counts.get(row["confidence"], 0) + 1
        if row["H1_lineage"] == "filtered_out":
            if row["confidence"] == "UNCLASSIFIED":
                filtered_unclassifiable += 1
            else:
                filtered_pre_pandemic += 1

    summary_path = out_dir / "classification_summary.txt"
    with summary_path.open("w") as handle:
        handle.write("=== H1N1 Lineage Classification Summary ===\n\n")
        handle.write(f"Total sequences: {len(rows)}\n")
        handle.write(f"  pdm09: {counts['pdm09']}\n")
        handle.write(f"  seasonal: {counts['seasonal']}\n")
        handle.write(f"  filtered_out: {counts['filtered_out']}\n")
        handle.write(f"    pre-{args.pdm_min_year} pdm09: {filtered_pre_pandemic}\n")
        handle.write(f"    unclassifiable (divergence gate): {filtered_unclassifiable}\n\n")
        handle.write("Confidence distribution:\n")
        for key in ["HIGH", "MODERATE", "AMBIGUOUS", "LOW", "UNCLASSIFIED"]:
            handle.write(f"  {key}: {conf_counts[key]} sequences\n")
        handle.write("\nThresholds:\n")
        handle.write(f"  max per-site divergence: {args.max_divergence_frac}\n")
        handle.write(f"  max absolute divergence (span unknown): {args.max_divergence_abs}\n")
        handle.write(f"  pdm09 minimum year: {args.pdm_min_year}\n\n")
        handle.write("References:\n")
        handle.write(f"  pdm09: {args.pdm_tsv}\n")
        handle.write(f"  seasonal: {args.seasonal_tsv}\n")

    print(
        f"[classify_h1n1_lineage] pdm09={counts['pdm09']} seasonal={counts['seasonal']} "
        f"filtered_out={counts['filtered_out']} "
        f"(HIGH={conf_counts['HIGH']} MODERATE={conf_counts['MODERATE']} "
        f"AMBIGUOUS={conf_counts['AMBIGUOUS']} LOW={conf_counts['LOW']} "
        f"UNCLASSIFIED={conf_counts['UNCLASSIFIED']})"
    )


def write_filter_out(path: Path, filtered_rows: list[dict[str, str]], metadata_path: str | None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lineage_frame = pd.DataFrame(filtered_rows, columns=LINEAGE_FIELDS)
    if metadata_path and Path(metadata_path).exists() and not lineage_frame.empty:
        meta = pd.read_csv(metadata_path, dtype=str, keep_default_na=False)
        if "Isolate_Id" in meta.columns:
            # Full metadata for each filtered isolate, annotated with the reason.
            extra = lineage_frame[
                [
                    "Isolate_Id",
                    "pdm_divergence",
                    "seasonal_divergence",
                    "pdm_divergence_frac",
                    "seasonal_divergence_frac",
                    "filter_reason",
                ]
            ]
            merged = meta.merge(extra, on="Isolate_Id", how="inner")
            merged.to_csv(path, index=False)
            return
    lineage_frame.to_csv(path, index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Classify H1N1 isolates as seasonal or pdm09 from HA alignment divergence.")
    parser.add_argument("--pdm-tsv", required=True, help="Nextclade TSV of HA aligned to the pdm09 reference")
    parser.add_argument("--seasonal-tsv", required=True, help="Nextclade TSV of HA aligned to the seasonal reference")
    parser.add_argument("--out", required=True, help="Output lineage assignment CSV")
    parser.add_argument("--filter-out", default="lineage_filter_out.csv", help="CSV of filtered-out (pre-pandemic pdm09) isolates")
    parser.add_argument("--metadata", default="", help="Merged metadata CSV (for Year and the filtered-out report)")
    parser.add_argument("--max-divergence-frac", type=float, default=0.15, help="Max per-site divergence for a real lineage match")
    parser.add_argument("--max-divergence-abs", type=float, default=250.0, help="Max absolute divergence when aligned length is unknown")
    parser.add_argument("--pdm-min-year", type=int, default=2009, help="Filter out pdm09-classified isolates collected before this year")
    return parser.parse_args()


if __name__ == "__main__":
    main()
