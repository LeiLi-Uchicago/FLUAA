from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


DATE_RE = re.compile(r"^(\d{4})(?:[-/](\d{1,2})(?:[-/](\d{1,2}))?)?$")
NA_SUBTYPE_RE = re.compile(r"N([1-9])\b", re.IGNORECASE)
STRAIN_YEAR_RE = re.compile(r"/(\d{4})$")


def normalize_date(value: object) -> tuple[str, str]:
    if value is None or pd.isna(value):
        return "", ""
    if isinstance(value, pd.Timestamp):
        return f"{value.year:04d}", f"{value.month:02d}"
    text = str(value).strip()
    if not text or text.upper() in {"NA", "NAN", "NONE", "NULL"}:
        return "", ""
    match = DATE_RE.match(text)
    if match:
        year, month, _day = match.groups()
        return year, f"{int(month):02d}" if month else ""
    parsed = pd.to_datetime(text, errors="coerce")
    if pd.isna(parsed):
        return "", ""
    return f"{parsed.year:04d}", f"{parsed.month:02d}"


def normalize_date_for_strain(value: object, strain: object) -> tuple[str, str]:
    year, month = normalize_date(value)
    expected_year = strain_year(strain)
    if expected_year and year and expected_year != year:
        return expected_year, ""
    if expected_year and not year:
        return expected_year, ""
    return year, month


def strain_year(strain: object) -> str:
    if strain is None or pd.isna(strain):
        return ""
    match = STRAIN_YEAR_RE.search(str(strain).strip())
    return match.group(1) if match else ""


def read_metadata_files(input_dir: Path, id_column: str, date_column: str) -> pd.DataFrame:
    files = sorted(
        path
        for pattern in ("*.xls", "*.xlsx", "*.csv", "*.tsv")
        for path in input_dir.glob(pattern)
        if not is_hidden_or_sidecar(path)
    )
    if not files:
        raise FileNotFoundError(f"No .xls/.xlsx/.csv/.tsv metadata files found in {input_dir}")

    frames: list[pd.DataFrame] = []
    for order, path in enumerate(files):
        if path.suffix.lower() == ".csv":
            frame = pd.read_csv(path, dtype=object)
        elif path.suffix.lower() == ".tsv":
            frame = pd.read_csv(path, sep="\t", dtype=object)
        else:
            try:
                frame = pd.read_excel(path, sheet_name=0, dtype=object, engine=excel_engine_for(path))
            except ImportError as exc:
                raise RuntimeError(
                    f"Reading {path.name} requires an Excel engine. Install xlrd for .xls files "
                    "or use the provided conda environment."
                ) from exc
            except ValueError as exc:
                raise RuntimeError(
                    f"Could not read metadata file {path}. If this is a hidden AppleDouble file like '._*.xls', "
                    "it should be ignored by the pipeline; otherwise check that the file is a valid Excel workbook."
                ) from exc
        frame["_source_file"] = path.name
        frame["_source_order"] = order
        frame["_row_order"] = range(len(frame))
        frames.append(frame)

    data = pd.concat(frames, ignore_index=True)
    for column in [id_column, "Isolate_Name", date_column, "Subtype", "Lineage", "Update_Date", "Submission_Date"]:
        if column not in data.columns:
            data[column] = ""

    data[id_column] = data[id_column].map(_clean_text)
    data["Isolate_Name"] = data["Isolate_Name"].map(_clean_text)
    data["strain"] = data["Isolate_Name"]
    years_months = [
        normalize_date_for_strain(date_value, strain_value)
        for date_value, strain_value in zip(data[date_column], data["strain"])
    ]
    data["Year"] = [item[0] for item in years_months]
    data["Month"] = [item[1] for item in years_months]
    data["NA_subtype"] = data["Subtype"].map(parse_na_subtype)

    data["_update_sort"] = pd.to_datetime(data["Update_Date"], errors="coerce")
    data["_submission_sort"] = pd.to_datetime(data["Submission_Date"], errors="coerce")
    data = data.sort_values(
        [id_column, "_update_sort", "_submission_sort", "_source_order", "_row_order"],
        ascending=[True, False, False, True, True],
        na_position="last",
    )
    data = data.drop_duplicates(subset=[id_column], keep="first")
    data = data.rename(columns={"_source_file": "source_file"})
    return data.drop(columns=["_update_sort", "_submission_sort", "_source_order", "_row_order"])


def parse_na_subtype(value: object) -> str:
    if value is None or pd.isna(value):
        return "unknown"
    match = NA_SUBTYPE_RE.search(str(value))
    return f"N{match.group(1)}" if match else "unknown"


def _clean_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def is_hidden_or_sidecar(path: Path) -> bool:
    name = path.name
    return name.startswith(".") or name.startswith("._")


def excel_engine_for(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".xls":
        return "xlrd"
    if suffix == ".xlsx":
        return "openpyxl"
    raise ValueError(f"Unsupported Excel extension: {path}")
