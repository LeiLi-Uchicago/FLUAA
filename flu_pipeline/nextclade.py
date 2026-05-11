from __future__ import annotations

import json
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pandas as pd


GENE_NAMES = {"HA", "NA", "MP", "NP", "NS", "PA", "PB1", "PB2"}

HA_COLUMNS_BY_SUBTYPE = {
    "H1N1": ["HA_clade", "HA_proposedSubclade", "HA_subclade", "HA_short_clade", "HA_legacy_clade"],
    "H3N2": ["HA_clade", "HA_proposedSubclade", "HA_subclade", "HA_short_clade", "HA_legacy_clade"],
    "B_VIC": ["HA_clade", "HA_proposedSubclade", "HA_subclade", "HA_legacy_clade"],
    "B_YAM": ["HA_clade", "HA_legacy_clade_yam"],
}

NEXTCLADE_COLUMN_ALIASES = {
    "clade": ["clade", "nextcladePangoLineage", "Nextclade_pango"],
    "proposedSubclade": ["proposedSubclade", "proposed_subclade", "proposedSubcladeName"],
    "subclade": ["subclade", "sub_clade"],
    "short_clade": ["short-clade", "short_clade", "shortClade", "short_clade_label"],
    "legacy_clade": ["legacy-clade", "legacy_clade", "legacyClade", "legacyCladeName"],
    "legacy_clade_yam": ["legacy-clade", "legacy_clade_yam", "legacyCladeYam", "legacy_clade", "legacyClade"],
}


def read_nextclade_tsv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def seq_name_column(frame: pd.DataFrame) -> str:
    for column in ["seqName", "seq_name", "name", "sequenceName"]:
        if column in frame.columns:
            return column
    raise ValueError(f"Cannot find sequence-name column in Nextclade TSV columns: {list(frame.columns)}")


def resolve_column(frame: pd.DataFrame, logical_name: str) -> str | None:
    for candidate in NEXTCLADE_COLUMN_ALIASES.get(logical_name, [logical_name]):
        if candidate in frame.columns:
            return candidate
    return None


def iter_ndjson(path: Path) -> Iterator[dict[str, Any]]:
    if not path.exists():
        return
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line:
                yield json.loads(line)


def find_nested_lists(obj: Any, key_contains: str) -> Iterator[list[Any]]:
    needle = key_contains.lower()
    if isinstance(obj, dict):
        for key, value in obj.items():
            if needle in key.lower() and isinstance(value, list):
                yield value
            yield from find_nested_lists(value, key_contains)
    elif isinstance(obj, list):
        for item in obj:
            yield from find_nested_lists(item, key_contains)


def object_seq_name(obj: dict[str, Any]) -> str:
    for key in ["seqName", "seq_name", "name", "sequenceName"]:
        value = obj.get(key)
        if isinstance(value, str):
            return value
    return ""


def gene_from_nextclade_dir_name(name: str) -> str:
    text = name.upper()
    if text in GENE_NAMES:
        return text
    for gene in sorted(GENE_NAMES, key=len, reverse=True):
        if text.endswith(f"_{gene}"):
            return gene
    return text
