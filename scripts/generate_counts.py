#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Iterable

import pandas as pd
from Bio.Align import PairwiseAligner

from flu_pipeline.fasta import parse_isolate_from_seq_name
from flu_pipeline.nextclade import gene_from_nextclade_dir_name, iter_ndjson, object_seq_name


OBSERVED_AA = set("ACDEFGHIKLMNPQRSTVWY*")
AMBIGUOUS_AA = set("XBZJ")
MISSING_AA = set("?.")
VALID_CODON = set("ACGTU")
VALID_NT = set("ACGTRYMKSWBDHVN")
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

def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "insertions").mkdir(exist_ok=True)

    metadata = pd.read_csv(args.metadata, dtype=str, keep_default_na=False)
    metadata_by_id = metadata.set_index("Isolate_Id", drop=False).to_dict(orient="index")
    clade_columns = [col for col in args.clade_columns.split(",") if col and col in metadata.columns]

    events = collect_insertion_events(args.nextclade_dirs, metadata_by_id)
    insertion_support = summarize_insertions(events)
    write_insertion_outputs(outdir, events, insertion_support)
    write_supported_insertion_summary(outdir, insertion_support, args.insertion_min_support)

    counts_by_protein: dict[str, Counter[tuple[str, ...]]] = {}
    for directory in map(Path, args.nextclade_dirs):
        aligned_nt_sequences = read_fasta_sequences(directory / "aligned.fasta")
        for protein, records in iter_translation_record_groups(directory / "translations"):
            start = time.monotonic()
            print(f"[generate_counts] counting {protein} from {directory.name}", file=sys.stderr, flush=True)
            codon_lookup = build_codon_lookup(records, aligned_nt_sequences)
            counts = aggregate_translation(
                protein,
                records.items(),
                codon_lookup,
                metadata_by_id,
                clade_columns,
            )
            counts_by_protein.setdefault(protein, Counter()).update(counts)
            elapsed = time.monotonic() - start
            print(f"[generate_counts] counted {protein} from {directory.name} in {elapsed:.1f}s", file=sys.stderr, flush=True)

    for protein, counts in sorted(counts_by_protein.items()):
        write_count_tables_for_protein(protein, counts, outdir, clade_columns)
        print(f"[generate_counts] wrote merged {protein}", file=sys.stderr, flush=True)

    remove_empty_sqlite_from_previous_runs(outdir)


def aggregate_translation(
    protein: str,
    records: Iterable[tuple[str, str]],
    codon_lookup: dict[str, list[tuple[str, str]]],
    metadata_by_id: dict[str, dict[str, str]],
    clade_columns: list[str],
) -> Counter[tuple[str, ...]]:
    counts: Counter[tuple[str, ...]] = Counter()
    for seq_name, aa_sequence in records:
        isolate_id = parse_isolate_from_seq_name(seq_name)
        meta = metadata_by_id.get(isolate_id, {})
        codons = codon_lookup.get(seq_name, [])
        for index, aa in enumerate(aa_sequence, start=1):
            codon, codon_source = codons[index - 1] if index - 1 < len(codons) else ("NA", "unmapped")
            amino_acid, codon, codon_status = normalize_amino_acid_and_codon(
                aa,
                codon,
            )
            add_observation(counts, protein, str(index), amino_acid, codon, codon_status, codon_source, meta, clade_columns)
    return counts


def add_observation(
    batch: Counter[tuple[str, ...]],
    protein: str,
    position: str,
    amino_acid: str,
    codon: str,
    codon_status: str,
    codon_source: str,
    meta: dict[str, str],
    clade_columns: list[str],
) -> None:
    year = unknown_if_empty(meta.get("Year", ""))
    month = unknown_if_empty(meta.get("Month", ""))
    rows = [("", "")]
    for column in clade_columns:
        rows.append((column, unknown_if_empty(meta.get(column, ""))))
    for group_column, group_value in rows:
        batch[(protein, position, amino_acid, codon, codon_status, codon_source, year, month, group_column, group_value)] += 1


def classify_aa(aa: str, codon: str) -> tuple[str, str, str]:
    aa = aa.upper()
    if aa == "-":
        return "-", "deletion", "DEL"
    if aa in MISSING_AA:
        return aa, "missing", "NA"
    if aa in AMBIGUOUS_AA:
        return aa, "ambiguous", "NA"
    if aa not in OBSERVED_AA:
        return aa, "ambiguous", "NA"
    codon = codon.upper().replace("U", "T")
    if len(codon) != 3 or any(base not in VALID_CODON for base in codon):
        codon = "NA"
    return aa, "observed", codon


def normalize_amino_acid_and_codon(aa: str, codon: str) -> tuple[str, str, str]:
    aa = aa.upper()
    clean_codon, codon_status = normalize_codon_with_status(codon)
    if aa == "-":
        return "-", "DEL", "deletion"
    if aa in MISSING_AA:
        return "Unknown", "NA", "missing_aa"
    if not aa:
        return "Unknown", "NA", "missing_aa"
    if aa in AMBIGUOUS_AA:
        if codon_status in {"valid_codon", "codon_ambiguous"}:
            return aa, clean_codon, "ambiguous_aa"
        return aa, "NA", codon_status
    if aa not in OBSERVED_AA:
        return aa, "NA", "ambiguous_aa"
    if codon_status != "valid_codon":
        return aa, "NA", codon_status
    if CODON_TABLE.get(clean_codon) != aa:
        return aa, clean_codon, "codon_aa_mismatch"
    return aa, clean_codon, "observed_exact"


def normalize_codon(codon: str) -> str:
    return normalize_codon_with_status(codon)[0]


def normalize_codon_with_status(codon: str) -> tuple[str, str]:
    codon = codon.upper().replace("U", "T")
    if codon in {"", "NA", "NAN", "NONE", "NULL"}:
        return "NA", "codon_unavailable"
    if codon == "DEL":
        return "DEL", "deletion"
    if len(codon) != 3:
        return "NA", "codon_incomplete"
    if any(base not in VALID_NT for base in codon):
        if len(codon) == 3:
            return "NA", "codon_ambiguous"
    if any(base not in {"A", "C", "G", "T"} for base in codon):
        return codon, "codon_ambiguous"
    return codon, "valid_codon"


def unknown_if_empty(value: object) -> str:
    text = "" if value is None else str(value).strip()
    if not text or text.upper() in {"NA", "NAN", "NONE", "NULL"}:
        return "Unknown"
    return text


def codons_from_aligned_ranges(nt_sequence: str, ranges: list[tuple[int, int]], frame_offset: int = 0) -> list[str]:
    ungapped_sequence = "".join(base for base in nt_sequence.upper() if base not in {"-", ".", " "})
    segment = "".join(ungapped_sequence[begin:end] for begin, end in ranges)[frame_offset:]
    return [segment[i : i + 3] for i in range(0, len(segment) - 2, 3)]


def codons_for_aligned_peptide(aa_sequence: str, raw_codons: list[str]) -> list[str]:
    codons: list[str] = []
    raw_index = 0
    for aa in aa_sequence:
        if aa == "-":
            codons.append("DEL")
            continue
        if raw_index < len(raw_codons):
            codons.append(raw_codons[raw_index])
            raw_index += 1
        else:
            codons.append("NA")
    return codons


def best_codons_for_aligned_peptide(aa_sequence: str, candidate_codons: list[list[str]]) -> list[str]:
    if not candidate_codons:
        return []
    best = min(
        (codons_for_aligned_peptide(aa_sequence, codons) for codons in candidate_codons),
        key=lambda codons: codon_match_penalty(aa_sequence, codons),
    )
    return best


def codon_match_penalty(aa_sequence: str, codons: list[str]) -> int:
    penalty = 0
    for index, aa in enumerate(aa_sequence):
        aa = aa.upper()
        if aa == "-" or aa in MISSING_AA or aa in AMBIGUOUS_AA or aa not in OBSERVED_AA:
            continue
        codon = normalize_codon(codons[index] if index < len(codons) else "NA")
        if codon == "NA":
            penalty += 10
        elif CODON_TABLE.get(codon) != aa:
            penalty += 1000
    return penalty


def read_fasta_sequences(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    return dict(iter_simple_fasta(path))


def iter_translation_record_groups(translations_dir: Path) -> Iterable[tuple[str, dict[str, str]]]:
    files = {path.stem: path for path in sorted(translations_dir.glob("*.fasta"))}
    for ignored in ["PA-X"]:
        files.pop(ignored, None)

    if {"SigPep", "HA1", "HA2"}.issubset(files):
        component_records = {
            "SigPep": read_fasta_sequences(files["SigPep"]),
            "HA1": read_fasta_sequences(files["HA1"]),
            "HA2": read_fasta_sequences(files["HA2"]),
        }
        joined_records = joined_ha_records(component_records)
        yield "HA", joined_records
        for key in ["SigPep", "HA1", "HA2"]:
            files.pop(key, None)

    for protein, path in sorted(files.items()):
        yield protein, read_fasta_sequences(path)


def joined_ha_records(component_records: dict[str, dict[str, str]]) -> dict[str, str]:
    sigpep = component_records["SigPep"]
    ha1 = component_records["HA1"]
    ha2 = component_records["HA2"]
    names = sorted(set(sigpep) | set(ha1) | set(ha2))
    return {name: sigpep.get(name, "") + ha1.get(name, "") + ha2.get(name, "") for name in names}


def build_codon_lookup(
    records: dict[str, str],
    aligned_nt_sequences: dict[str, str],
) -> dict[str, list[tuple[str, str]]]:
    aligned_by_id = {
        parse_isolate_from_seq_name(seq_name): sequence
        for seq_name, sequence in aligned_nt_sequences.items()
    }
    lookup: dict[str, list[tuple[str, str]]] = {}
    for seq_name, aa_sequence in records.items():
        isolate_id = parse_isolate_from_seq_name(seq_name)
        aligned_nt = aligned_nt_sequences.get(seq_name) or aligned_by_id.get(isolate_id)
        lookup[seq_name] = choose_codon_mapping(aa_sequence, aligned_nt)
    return lookup


def choose_codon_mapping(aa_sequence: str, aligned_nt: str | None) -> list[tuple[str, str]]:
    if not aligned_nt:
        return [(("NA" if aa != "-" else "DEL"), "unmapped") for aa in aa_sequence]
    codons = map_codons_from_nt_to_gapped_aa(aligned_nt, aa_sequence)
    return [(codon, "deletion" if codon == "DEL" else "nextclade_aligned_nt") for codon in codons]


def map_codons_from_nt_to_gapped_aa(nt_sequence: str, aa_sequence: str) -> list[str]:
    aa_sequence = aa_sequence.upper()
    observed_aa = "".join(aa for aa in aa_sequence if aa != "-")
    if not observed_aa:
        return ["DEL" for _aa in aa_sequence]

    candidates: list[list[str]] = []
    for _strand, frame, codons, translated in translated_frame_candidates(nt_sequence):
        exact = find_observed_aa_start(translated, observed_aa)
        if exact is not None:
            candidates.append(gap_codons_for_aa_sequence(aa_sequence, codons, exact))

    if not candidates:
        aligner = protein_aligner()
        for _strand, frame, codons, translated in translated_frame_candidates(nt_sequence):
            greedy = greedy_observed_aa_mapping(observed_aa, translated)
            if greedy:
                candidates.append(gap_codons_from_query_mapping(aa_sequence, codons, greedy))
            mapped = align_observed_aa_to_translated_frame(aligner, observed_aa, translated)
            if mapped:
                candidates.append(gap_codons_from_query_mapping(aa_sequence, codons, mapped))

    if not candidates:
        return ["NA" if aa != "-" else "DEL" for aa in aa_sequence]
    return min(candidates, key=lambda codons: codon_match_penalty(aa_sequence, codons))


def translated_frame_candidates(nt_sequence: str) -> Iterable[tuple[str, int, list[str], str]]:
    clean_nt = clean_nt_sequence(nt_sequence)
    for strand, sequence in [("+", clean_nt), ("-", reverse_complement(clean_nt))]:
        for frame in range(3):
            codons = [sequence[index : index + 3] for index in range(frame, len(sequence) - 2, 3)]
            yield strand, frame, codons, translate_codons(codons)


def clean_nt_sequence(nt_sequence: str) -> str:
    return "".join(base for base in nt_sequence.upper().replace("U", "T") if base.isalpha())


def translate_codons(codons: list[str]) -> str:
    return "".join(translate_codon(codon) for codon in codons)


def translate_codon(codon: str) -> str:
    codon, status = normalize_codon_with_status(codon)
    if status == "valid_codon":
        return CODON_TABLE.get(codon, "X")
    if status == "codon_ambiguous":
        return "X"
    return "X"


def find_observed_aa_start(translated: str, observed_aa: str) -> int | None:
    if len(observed_aa) > len(translated):
        return None
    best_start = None
    best_score = -1
    for start in range(0, len(translated) - len(observed_aa) + 1):
        score = 0
        ok = True
        for query_aa, translated_aa in zip(observed_aa, translated[start : start + len(observed_aa)]):
            if aa_translation_match(query_aa, translated_aa):
                score += 1
            else:
                ok = False
                break
        if ok and score > best_score:
            best_start = start
            best_score = score
    return best_start


def aa_translation_match(query_aa: str, translated_aa: str) -> bool:
    return query_aa == translated_aa or query_aa == "X" or translated_aa == "X"


def gap_codons_for_aa_sequence(aa_sequence: str, codons: list[str], start: int) -> list[str]:
    output: list[str] = []
    codon_index = start
    for aa in aa_sequence:
        if aa == "-":
            output.append("DEL")
            continue
        output.append(codons[codon_index] if codon_index < len(codons) else "NA")
        codon_index += 1
    return output


def protein_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.01
    aligner.end_gap_score = 0
    return aligner


def align_observed_aa_to_translated_frame(
    aligner: PairwiseAligner,
    observed_aa: str,
    translated: str,
) -> dict[int, int]:
    if not observed_aa or not translated:
        return {}
    alignment = aligner.align(translated, observed_aa)[0]
    target_blocks, query_blocks = alignment.aligned
    mapping: dict[int, int] = {}
    for (target_start, target_end), (query_start, query_end) in zip(target_blocks, query_blocks):
        block_size = min(target_end - target_start, query_end - query_start)
        for offset in range(block_size):
            mapping[query_start + offset] = target_start + offset
    return mapping


def greedy_observed_aa_mapping(observed_aa: str, translated: str) -> dict[int, int]:
    mapping: dict[int, int] = {}
    target_index = 0
    for query_index, query_aa in enumerate(observed_aa):
        found = None
        for index in range(target_index, len(translated)):
            if aa_translation_match(query_aa, translated[index]):
                found = index
                break
        if found is None:
            if target_index >= len(translated):
                return {}
            found = target_index
        mapping[query_index] = found
        target_index = found + 1
    return mapping


def gap_codons_from_query_mapping(
    aa_sequence: str,
    codons: list[str],
    query_to_target: dict[int, int],
) -> list[str]:
    output: list[str] = []
    observed_index = 0
    for aa in aa_sequence:
        if aa == "-":
            output.append("DEL")
            continue
        target_index = query_to_target.get(observed_index)
        output.append(codons[target_index] if target_index is not None and target_index < len(codons) else "NA")
        observed_index += 1
    return output


def cds_candidate_ranges(obj: dict) -> dict[str, list[list[tuple[int, int]]]]:
    ranges: dict[str, list[list[tuple[int, int]]]] = {}
    annotation = obj.get("annotation")
    if not isinstance(annotation, dict):
        return ranges
    genes = annotation.get("genes")
    if not isinstance(genes, list):
        return ranges
    for gene in genes:
        if not isinstance(gene, dict):
            continue
        name = str(gene.get("name") or gene.get("id") or "")
        if not name:
            continue
        cds_segments = []
        for cds in gene.get("cdses") or []:
            if not isinstance(cds, dict):
                continue
            cds_segments.extend(cds.get("segments") or [])
        candidates = []
        for key in ["range", "rangeLocal"]:
            segment_ranges = ranges_from_cds_segments(cds_segments, key)
            if segment_ranges and segment_ranges not in candidates:
                candidates.append(segment_ranges)
        if candidates:
            ranges[name] = candidates
    return ranges


def ranges_from_cds_segments(segments: list, range_key: str) -> list[tuple[int, int]]:
    segment_ranges: list[tuple[int, int]] = []
    for segment in segments:
        if not isinstance(segment, dict):
            continue
        range_obj = segment.get(range_key)
        if not isinstance(range_obj, dict):
            continue
        begin = range_obj.get("begin")
        end = range_obj.get("end")
        if isinstance(begin, int) and isinstance(end, int) and end > begin:
            segment_ranges.append((begin, end))
    return sorted(segment_ranges)


def reverse_complement(sequence: str) -> str:
    complement = str.maketrans(
        "ACGTUacgtuNnRYMKSWBDHVrymkswbdhv-",
        "TGCAAtgcaaNnYRKMWSVHDByrkmwsvhdb-",
    )
    return sequence.translate(complement)[::-1]


def single_cds_codon_lookup(
    cds_codons: dict[str, dict[str, list[list[str]]]],
    protein: str,
    records: dict[str, str],
) -> dict[str, list[str]]:
    return {
        seq_name: best_codons_for_aligned_peptide(records.get(seq_name, ""), per_cds.get(protein, []))
        for seq_name, per_cds in cds_codons.items()
        if seq_name in records
    }


def joined_codon_lookup(
    cds_codons: dict[str, dict[str, list[list[str]]]],
    component_records: dict[str, dict[str, str]],
) -> dict[str, list[str]]:
    joined: dict[str, list[str]] = {}
    proteins = ["SigPep", "HA1", "HA2"]
    names = sorted(set().union(*(records.keys() for records in component_records.values())))
    for seq_name, per_cds in cds_codons.items():
        if seq_name not in names:
            continue
        codons: list[str] = []
        for protein in proteins:
            if seq_name in component_records[protein]:
                codons.extend(best_codons_for_aligned_peptide(component_records[protein][seq_name], per_cds.get(protein, [])))
        joined[seq_name] = codons
    return joined


def iter_simple_fasta(path: Path) -> Iterable[tuple[str, str]]:
    header = None
    chunks: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def collect_insertion_events(nextclade_dirs: list[str], metadata_by_id: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    events: list[dict[str, str]] = []
    for directory in map(Path, nextclade_dirs):
        fallback_gene = gene_from_nextclade_dir_name(directory.name)
        for obj in iter_ndjson(directory / "nextclade.ndjson"):
            seq_name = object_seq_name(obj)
            isolate_id = parse_isolate_from_seq_name(seq_name)
            meta = metadata_by_id.get(isolate_id, {})
            for item in insertion_items(obj):
                protein = str(item.get("cdsName") or item.get("cds") or item.get("gene") or fallback_gene).upper()
                position = str(item.get("pos") or item.get("position") or item.get("refPos") or item.get("left") or "")
                aa = str(item.get("aa") or item.get("query") or item.get("ins") or item.get("inserted") or "INS")
                codon = str(item.get("codon") or item.get("nuc") or item.get("nucSequence") or "NA").upper()
                if not position:
                    continue
                if len(codon) != 3 or any(base not in VALID_CODON for base in codon.replace("U", "T")):
                    codon = "NA"
                events.append(
                    {
                        "seq_name": seq_name,
                        "Isolate_Id": isolate_id,
                        "protein": protein,
                        "position_label": f"{position}ins",
                        "aa_state": aa,
                        "codon": codon.replace("U", "T"),
                        "Year": meta.get("Year", ""),
                        "Month": meta.get("Month", ""),
                    }
                )
    return events


def insertion_items(obj: dict) -> Iterable[dict]:
    for key in ["insertions", "aaInsertions"]:
        value = obj.get(key)
        if isinstance(value, list):
            for item in value:
                if isinstance(item, dict):
                    yield item


def summarize_insertions(events: list[dict[str, str]]) -> Counter[tuple[str, str, str, str]]:
    support: Counter[tuple[str, str, str, str]] = Counter()
    seen: set[tuple[str, str, str, str, str]] = set()
    for event in events:
        key = (event["protein"], event["position_label"], event["aa_state"], event["codon"])
        strain_key = (*key, event["Isolate_Id"])
        if strain_key not in seen:
            support[key] += 1
            seen.add(strain_key)
    return support


def write_insertion_outputs(outdir: Path, events: list[dict[str, str]], support: Counter[tuple[str, str, str, str]]) -> None:
    event_path = outdir / "insertions" / "insertion_events.csv"
    event_fields = ["seq_name", "Isolate_Id", "protein", "position_label", "aa_state", "codon", "Year", "Month"]
    with event_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=event_fields)
        writer.writeheader()
        writer.writerows(events)

    summary_path = outdir / "insertions" / "insertion_summary.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["protein", "position_label", "aa_state", "codon", "support"])
        writer.writeheader()
        for (protein, position_label, aa_state, codon), count in sorted(support.items()):
            writer.writerow(
                {
                    "protein": protein,
                    "position_label": position_label,
                    "aa_state": aa_state,
                    "codon": codon,
                    "support": count,
                }
            )


def write_supported_insertion_summary(
    outdir: Path,
    support: Counter[tuple[str, str, str, str]],
    min_support: int,
) -> None:
    path = outdir / "insertions" / "supported_insertions.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["protein", "position_label", "aa_state", "codon", "support"])
        writer.writeheader()
        for (protein, position_label, aa_state, codon), count in sorted(support.items()):
            if count >= min_support:
                writer.writerow(
                    {
                        "protein": protein,
                        "position_label": position_label,
                        "aa_state": aa_state,
                        "codon": codon,
                        "support": count,
                    }
                )


def write_count_tables_for_protein(
    protein: str,
    counts: Counter[tuple[str, ...]],
    outdir: Path,
    clade_columns: list[str],
) -> None:
    protein_dir = outdir / protein
    protein_dir.mkdir(parents=True, exist_ok=True)
    write_year_month_table(counts, protein_dir / "aa_usage_by_Year_Month.csv")
    for column in clade_columns:
        write_group_table(counts, column, protein_dir / f"aa_usage_by_{column}.csv")


def write_year_month_table(counts: Counter[tuple[str, ...]], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Protein", "Position", "Year", "Month", "AminoAcid", "Codon", "CodonStatus", "CodonSource", "Count"])
        writer.writeheader()
        rows = [
            (protein, position, year, month, amino_acid, codon, codon_status, codon_source, count)
            for (protein, position, amino_acid, codon, codon_status, codon_source, year, month, group_column, _group_value), count in counts.items()
            if group_column == ""
        ]
        for row in sorted(rows, key=count_sort_key):
            writer.writerow(
                {
                    "Protein": row[0],
                    "Position": row[1],
                    "Year": row[2],
                    "Month": row[3],
                    "AminoAcid": row[4],
                    "Codon": row[5],
                    "CodonStatus": row[6],
                    "CodonSource": row[7],
                    "Count": row[8],
                }
            )


def write_group_table(counts: Counter[tuple[str, ...]], column: str, path: Path) -> None:
    fieldnames = ["Protein", "Position", column, "Year", "Month", "AminoAcid", "Codon", "CodonStatus", "CodonSource", "Count"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        rows = [
            (protein, position, group_value, year, month, amino_acid, codon, codon_status, codon_source, count)
            for (protein, position, amino_acid, codon, codon_status, codon_source, year, month, group_column, group_value), count in counts.items()
            if group_column == column
        ]
        for row in sorted(rows, key=group_count_sort_key):
            writer.writerow(
                {
                    "Protein": row[0],
                    "Position": row[1],
                    column: row[2],
                    "Year": row[3],
                    "Month": row[4],
                    "AminoAcid": row[5],
                    "Codon": row[6],
                    "CodonStatus": row[7],
                    "CodonSource": row[8],
                    "Count": row[9],
                }
            )


def count_sort_key(row: tuple[object, ...]) -> tuple[int, str, str, str, str, str, str, str]:
    _protein, position, year, month, amino_acid, codon, codon_status, codon_source, _count = row
    return position_sort_value(str(position)), str(year), str(month), str(amino_acid), str(codon), str(codon_status), str(codon_source), str(position)


def group_count_sort_key(row: tuple[object, ...]) -> tuple[int, str, str, str, str, str, str, str, str]:
    _protein, position, group_value, year, month, amino_acid, codon, codon_status, codon_source, _count = row
    return position_sort_value(str(position)), str(group_value), str(year), str(month), str(amino_acid), str(codon), str(codon_status), str(codon_source), str(position)


def position_sort_value(position: str) -> int:
    try:
        return int(position)
    except ValueError:
        return 10**9


def remove_empty_sqlite_from_previous_runs(outdir: Path) -> None:
    path = outdir / "counts.sqlite"
    if path.exists():
        path.unlink()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate per-protein AA/codon count tables from Nextclade outputs.")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--nextclade-dirs", nargs="+", required=True)
    parser.add_argument("--gene-fastas", nargs="*", default=[])
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--insertion-min-support", type=int, default=2)
    parser.add_argument("--clade-columns", default="")
    return parser.parse_args()


if __name__ == "__main__":
    main()
