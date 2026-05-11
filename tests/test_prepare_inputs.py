from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

from flu_pipeline.fasta import parse_header
from flu_pipeline.metadata import normalize_date
from scripts.prepare_inputs import main as prepare_main


def test_parse_header() -> None:
    assert parse_header("A/Test/1/2020|EPI_ISL_1|HA") == ("A/Test/1/2020", "EPI_ISL_1", "HA")


def test_normalize_partial_dates() -> None:
    assert normalize_date("2020") == ("2020", "")
    assert normalize_date("2020-7") == ("2020", "07")
    assert normalize_date("2020-07-31") == ("2020", "07")
    assert normalize_date("not-a-date") == ("", "")


def test_prepare_inputs_merges_and_deduplicates(tmp_path: Path, monkeypatch) -> None:
    input_dir = tmp_path / "H1N1"
    input_dir.mkdir()
    pd.DataFrame(
        [
            {"Isolate_Id": "EPI_ISL_1", "Isolate_Name": "A/Test/1/2020", "Collection_Date": "2020-07-01", "Update_Date": "2020-08-01"},
            {"Isolate_Id": "EPI_ISL_2", "Isolate_Name": "A/Test/2/2020", "Collection_Date": "2020", "Update_Date": "2020-08-01"},
        ]
    ).to_csv(input_dir / "meta.csv", index=False)
    (input_dir / "sample-NT.fasta").write_text(
        ">A/Test/1/2020|EPI_ISL_1|HA\n"
        "ATGAAA\n"
        ">A/Test/1/2020|EPI_ISL_1|HA\n"
        "ATGAAN\n"
        ">A/Test/2/2020|EPI_ISL_2|NA\n"
        "ATGCCC\n"
    )
    outdir = tmp_path / "out"
    monkeypatch.setattr(
        "sys.argv",
        [
            "prepare_inputs.py",
            "--input-dir",
            str(input_dir),
            "--outdir",
            str(outdir),
            "--subtype",
            "H1N1",
        ],
    )
    prepare_main()

    metadata = pd.read_csv(outdir / "metadata" / "merged_metadata.csv", dtype=str, keep_default_na=False)
    assert list(metadata["Isolate_Id"]) == ["EPI_ISL_1", "EPI_ISL_2"]
    assert list(metadata["Year"]) == ["2020", "2020"]
    assert (outdir / "fasta_by_gene" / "HA.fasta").read_text().count(">") == 1
    assert (outdir / "fasta_by_gene" / "NA.fasta").read_text().count(">") == 1
    with (outdir / "reports" / "fasta_duplicate_conflicts.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["kept_header"] == "A/Test/1/2020|EPI_ISL_1|HA"
