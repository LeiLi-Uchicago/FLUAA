#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AMBIGUOUS_AA = {"X", "B", "Z", "J"}
MISSING_AA = {"UNKNOWN", "?", "."}
VALID_NT = set("ACGTRYMKSWBDHVN")
CORE_COLUMNS = {"Protein", "Position", "Year", "Month", "AminoAcid", "Codon", "CodonStatus", "CodonSource", "Count"}


def main() -> None:
    args = parse_args()
    count_root = Path(args.count_root)
    report_path = Path(args.out)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    summary: Counter[tuple[str, str]] = Counter()
    source_summary: Counter[tuple[str, str, str]] = Counter()
    for path in selected_count_files(count_root):
        protein = path.parent.name
        for row in validate_file(path, protein):
            rows.append(row)
            summary[(protein, str(row["issue"]))] += int(row["Count"])
            source_summary[(protein, str(row["issue"]), str(row["CodonSource"]))] += int(row["Count"])

    fields = [
        "Protein",
        "issue",
        "Position",
        "Grouping",
        "GroupingValue",
        "Year",
        "Month",
        "AminoAcid",
        "Codon",
        "CodonStatus",
        "CodonSource",
        "TranslatedAA",
        "Count",
        "SourceTable",
    ]
    with report_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    summary_path = report_path.with_name(report_path.stem + "_summary.csv")
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Protein", "issue", "Count"])
        writer.writeheader()
        for (protein, issue), count in sorted(summary.items()):
            writer.writerow({"Protein": protein, "issue": issue, "Count": count})

    source_summary_path = report_path.with_name(report_path.stem + "_source_summary.csv")
    with source_summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Protein", "issue", "CodonSource", "Count"])
        writer.writeheader()
        for (protein, issue, codon_source), count in sorted(source_summary.items()):
            writer.writerow(
                {
                    "Protein": protein,
                    "issue": issue,
                    "CodonSource": codon_source,
                    "Count": count,
                }
            )

    print(f"Wrote {report_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {source_summary_path}")


def selected_count_files(count_root: Path) -> list[Path]:
    paths: list[Path] = []
    for protein_dir in sorted(path for path in count_root.iterdir() if path.is_dir() and path.name != "insertions"):
        preferred = protein_dir / "aa_usage_by_HA_clade.csv"
        if preferred.exists():
            paths.append(preferred)
            continue
        fallback = protein_dir / "aa_usage_by_Year_Month.csv"
        if fallback.exists():
            paths.append(fallback)
            continue
        paths.extend(sorted(protein_dir.glob("aa_usage_by_*.csv"))[:1])
    return paths


def validate_file(path: Path, protein: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        grouping_columns = [column for column in reader.fieldnames or [] if column not in CORE_COLUMNS]
        grouping = ";".join(grouping_columns)
        for row in reader:
            aa = row.get("AminoAcid", "")
            codon = row.get("Codon", "")
            codon_status = row.get("CodonStatus", "")
            codon_source = row.get("CodonSource", "")
            count = int(row.get("Count", "0") or 0)
            issue, translated = classify_issue(aa, codon, codon_status)
            rows.append(
                {
                    "Protein": protein,
                    "issue": issue,
                    "Position": row.get("Position", ""),
                    "Grouping": grouping,
                    "GroupingValue": ";".join(row.get(column, "") for column in grouping_columns),
                    "Year": row.get("Year", ""),
                    "Month": row.get("Month", ""),
                    "AminoAcid": aa,
                    "Codon": codon,
                    "CodonStatus": codon_status or infer_codon_status(aa, codon),
                    "CodonSource": codon_source,
                    "TranslatedAA": translated,
                    "Count": count,
                    "SourceTable": path.name,
                }
            )
    return rows


def classify_issue(aa: str, codon: str, codon_status: str = "") -> tuple[str, str]:
    aa = aa.strip().upper()
    codon = codon.strip().upper().replace("U", "T")
    codon_status = codon_status.strip().lower()

    if aa == "-":
        if codon == "DEL" and codon_status in {"", "deletion"}:
            return "valid_deletion", ""
        return "deletion_codon_not_DEL", ""

    if aa in MISSING_AA:
        if codon == "NA" and codon_status in {"", "missing_aa", "codon_unavailable"}:
            return "valid_missing_aa", ""
        return "missing_aa_codon_not_NA", ""

    if aa in AMBIGUOUS_AA:
        if codon == "NA" and codon_status in {"", "ambiguous_aa", "codon_unavailable"}:
            return "valid_ambiguous_aa", ""
        if is_codon_like(codon) and codon_status in {"ambiguous_aa", "codon_ambiguous"}:
            return "valid_ambiguous_aa", translate_if_standard(codon)
        return "ambiguous_aa_invalid_codon", ""

    if codon_status == "codon_unavailable":
        return "unrecovered_codon_for_observed_aa", ""
    if codon_status == "codon_ambiguous":
        return "ambiguous_codon_for_observed_aa", ""
    if codon_status == "codon_incomplete":
        return "incomplete_codon_for_observed_aa", ""
    if codon_status == "codon_aa_mismatch":
        translated = CODON_TABLE.get(codon, "")
        return "codon_mismatch", translated

    if codon == "NA":
        return "unrecovered_codon_for_observed_aa", ""
    if codon == "DEL":
        return "DEL_codon_for_observed_aa", ""
    if is_ambiguous_codon(codon):
        return "ambiguous_codon_for_observed_aa", ""
    translated = CODON_TABLE.get(codon)
    if translated is None:
        return "invalid_codon", ""
    if translated != aa:
        return "codon_mismatch", translated
    if codon_status and codon_status != "observed_exact":
        return f"status_conflict_{codon_status}", translated
    return "valid_observed_codon", translated


def is_codon_like(codon: str) -> bool:
    return len(codon) == 3 and all(base in VALID_NT for base in codon)


def is_ambiguous_codon(codon: str) -> bool:
    return is_codon_like(codon) and any(base not in {"A", "C", "G", "T"} for base in codon)


def translate_if_standard(codon: str) -> str:
    return CODON_TABLE.get(codon, "") if len(codon) == 3 else ""


def infer_codon_status(aa: str, codon: str) -> str:
    issue, _translated = classify_issue(aa, codon, "")
    if issue == "valid_observed_codon":
        return "observed_exact"
    if issue == "valid_deletion":
        return "deletion"
    if issue == "valid_missing_aa":
        return "missing_aa"
    if issue == "valid_ambiguous_aa":
        return "ambiguous_aa"
    return issue


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate amino-acid/codon consistency in count tables.")
    parser.add_argument("--count-root", required=True, help="Directory containing per-protein count folders, e.g. results/H1N1/count")
    parser.add_argument("--out", required=True, help="Path for detailed CSV report")
    return parser.parse_args()


if __name__ == "__main__":
    main()
