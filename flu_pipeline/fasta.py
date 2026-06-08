from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


AMBIGUOUS_NT = set("NRYKMSWBDHVU")
GAP_CHARS = set("-.")


@dataclass(frozen=True)
class FastaRecord:
    header: str
    strain: str
    isolate_id: str
    gene: str
    sequence: str


def parse_header(header: str) -> tuple[str, str, str]:
    text = header[1:] if header.startswith(">") else header
    parts = [part.strip() for part in text.strip().split("|")]
    if len(parts) == 3 and all(parts):
        return parts[0], parts[1], parts[2].upper()

    if len(parts) >= 3 and parts[0] and parts[-1]:
        for isolate_id in reversed(parts[1:-1]):
            if isolate_id:
                return parts[0], isolate_id, parts[-1].upper()

        raise ValueError(f"Expected FASTA header 'Isolate_Name|Isolate_Id|gene', got: {header!r}")
    raise ValueError(f"Expected FASTA header 'Isolate_Name|Isolate_Id|gene', got: {header!r}")


def iter_fasta(path: Path) -> Iterator[FastaRecord]:
    with path.open() as handle:
        header: str | None = None
        chunks: list[str] = []
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield _record_from_parts(header, chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield _record_from_parts(header, chunks)


def write_fasta_record(handle: TextIO, header: str, sequence: str, width: int = 80) -> None:
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), width):
        handle.write(sequence[start : start + width] + "\n")


def ungapped_len(sequence: str) -> int:
    return sum(1 for base in sequence.upper() if base not in GAP_CHARS)


def ambiguous_fraction(sequence: str) -> float:
    observed = [base for base in sequence.upper() if base not in GAP_CHARS]
    if not observed:
        return 1.0
    return sum(1 for base in observed if base not in {"A", "C", "G", "T"}) / len(observed)


def parse_isolate_from_seq_name(seq_name: str) -> str:
    parts = seq_name.split("|")
    if len(parts) >= 2 and parts[1].strip():
        return parts[1].strip()
    return seq_name.strip()


def _record_from_parts(header: str, chunks: list[str]) -> FastaRecord:
    strain, isolate_id, gene = parse_header(header)
    return FastaRecord(
        header=f"{strain}|{isolate_id}|{gene}",
        strain=strain,
        isolate_id=isolate_id,
        gene=gene,
        sequence="".join(chunks).upper(),
    )
