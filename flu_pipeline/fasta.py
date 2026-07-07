from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


AMBIGUOUS_NT = set("NRYKMSWBDHVU")
GAP_CHARS = set("-.")

# Canonical influenza gene/segment names, used to locate the gene field in a
# pipe-delimited FASTA header regardless of how many extra fields surround it.
GENE_NAMES = ("HA", "NA", "MP", "NP", "NS", "PA", "PB1", "PB2")


def normalize_gene(token: str) -> str | None:
    text = token.strip().upper()
    return text if text in GENE_NAMES else None


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
    if len(parts) < 3 or not parts[0]:
        raise ValueError(f"Expected FASTA header 'Isolate_Name|Isolate_Id|gene', got: {header!r}")

    strain = parts[0]

    # Locate the gene by matching a known gene name rather than assuming a fixed
    # position: GISAID exports may append a trailing segment id (e.g.
    # 'Name|Isolate_Id|HA|EPI351721') or carry extra fields before the gene.
    gene: str | None = None
    gene_index: int | None = None
    for index in range(len(parts) - 1, 0, -1):
        candidate = normalize_gene(parts[index])
        if candidate:
            gene, gene_index = candidate, index
            break

    if gene is None or gene_index is None:
        raise ValueError(f"Expected FASTA header 'Isolate_Name|Isolate_Id|gene', got: {header!r}")

    # The isolate id sits before the gene; prefer an explicit GISAID id.
    isolate_id = ""
    for part in parts[1:gene_index]:
        if part.upper().startswith("EPI_ISL"):
            isolate_id = part
            break
    if not isolate_id:
        for part in reversed(parts[1:gene_index]):
            if part:
                isolate_id = part
                break

    if not isolate_id:
        raise ValueError(f"Expected FASTA header 'Isolate_Name|Isolate_Id|gene', got: {header!r}")

    return strain, isolate_id, gene


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
