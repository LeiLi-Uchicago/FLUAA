from __future__ import annotations

from pathlib import Path

import pandas as pd

from scripts.annotate_metadata import main as annotate_main


def test_annotate_metadata_maps_ha_and_na_clades(tmp_path: Path, monkeypatch) -> None:
    metadata = tmp_path / "metadata.csv"
    pd.DataFrame([{"Isolate_Id": "EPI_ISL_1", "Isolate_Name": "A/Test/1/2020"}]).to_csv(metadata, index=False)
    ha_dir = tmp_path / "HA"
    na_dir = tmp_path / "NA"
    ha_dir.mkdir()
    na_dir.mkdir()
    pd.DataFrame(
        [
            {
                "seqName": "A/Test/1/2020|EPI_ISL_1|HA",
                "clade": "6B.1A",
                "proposedSubclade": "6B.1A.5a",
                "subclade": "sub",
                "short_clade": "short",
                "legacy_clade": "legacy",
            }
        ]
    ).to_csv(ha_dir / "nextclade.tsv", sep="\t", index=False)
    pd.DataFrame([{"seqName": "A/Test/1/2020|EPI_ISL_1|NA", "clade": "N1-a"}]).to_csv(
        na_dir / "nextclade.tsv", sep="\t", index=False
    )
    out = tmp_path / "annotated.csv"
    monkeypatch.setattr(
        "sys.argv",
        [
            "annotate_metadata.py",
            "--metadata",
            str(metadata),
            "--subtype",
            "H1N1",
            "--out",
            str(out),
            "--nextclade-dirs",
            str(ha_dir),
            str(na_dir),
        ],
    )
    annotate_main()
    result = pd.read_csv(out, dtype=str, keep_default_na=False)
    assert result.loc[0, "HA_clade"] == "6B.1A"
    assert result.loc[0, "HA_legacy_clade"] == "legacy"
    assert result.loc[0, "NA_clade"] == "N1-a"

