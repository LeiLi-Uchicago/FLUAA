#!/usr/bin/env python3
"""Standardize gap positions in HA translations to a canonical position.

For seasonal H1N1 HA, the stalk-region deletion (typically around position 147)
can align to different positions depending on local sequence similarity. This
script finds any gap cluster in a region around the canonical position and
moves it to the exact canonical position, ensuring all strains use the same
coordinate frame even though they carry the same underlying deletion.

Usage:
  standardize_gap_position.py --input translations.fasta --output normalized.fasta \\
    --protein HA --gap-position 147 --window 10
"""
from __future__ import annotations

import argparse
from pathlib import Path


def standardize_gaps(sequence: str, protein: str, gap_position: int, window: int) -> str:
    """Standardize gap positions in a sequence.

    For the specified protein, if there are gaps in the region [gap_position-window,
    gap_position+window], collect all gaps and reposition them to gap_position.
    Other proteins are left unchanged.
    """
    if protein != "HA":
        return sequence

    # Convert 1-indexed position to 0-indexed
    canonical_pos = gap_position - 1
    start = max(0, canonical_pos - window)
    end = min(len(sequence), canonical_pos + window + 1)

    # Extract region and check for gaps
    region = sequence[start:end]
    if "-" not in region:
        return sequence

    # Count gaps in the region
    gap_count = region.count("-")

    # Remove gaps from the region
    region_no_gaps = region.replace("-", "")

    # Rebuild: part before, canonical position with gaps, part after
    seq_before = sequence[:start]
    seq_after = sequence[end:]

    # Insert all gaps at the canonical position within the region
    region_with_gaps = (
        region_no_gaps[: canonical_pos - start]
        + "-" * gap_count
        + region_no_gaps[canonical_pos - start :]
    )

    return seq_before + region_with_gaps + seq_after


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)

    # Read FASTA
    lines = input_path.read_text().strip().split("\n")

    # Process sequences
    output_lines = []
    standardized_count = 0

    for i in range(0, len(lines), 2):
        if i + 1 >= len(lines):
            break

        header = lines[i]
        sequence = lines[i + 1]

        # Standardize gaps
        new_sequence = standardize_gaps(sequence, args.protein, args.gap_position, args.window)

        if new_sequence != sequence:
            standardized_count += 1

        output_lines.append(header)
        output_lines.append(new_sequence)

    # Write output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(output_lines) + "\n")

    print(
        f"[standardize_gap_position] {args.protein}: standardized {standardized_count} sequences "
        f"to gap position {args.gap_position}"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Standardize gap positions in protein translations (e.g., HA stalk deletions)."
    )
    parser.add_argument("--input", required=True, help="Input protein FASTA")
    parser.add_argument("--output", required=True, help="Output FASTA with standardized gaps")
    parser.add_argument(
        "--protein", default="HA", help="Protein to standardize (default: HA)"
    )
    parser.add_argument(
        "--gap-position",
        type=int,
        default=147,
        help="Canonical 1-indexed gap position (default: 147 for HA stalk)",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=10,
        help="Window around gap position to search (default: 10)",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
