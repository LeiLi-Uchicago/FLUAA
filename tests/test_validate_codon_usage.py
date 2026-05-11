from __future__ import annotations

from pathlib import Path

from scripts.validate_codon_usage import main as validate_main


def test_validate_codon_usage_writes_source_summary(tmp_path: Path, monkeypatch) -> None:
    count_root = tmp_path / "count"
    protein_dir = count_root / "HA"
    protein_dir.mkdir(parents=True)
    (protein_dir / "aa_usage_by_Year_Month.csv").write_text(
        "Protein,Position,Year,Month,AminoAcid,Codon,CodonStatus,CodonSource,Count\n"
        "HA,1,2020,01,M,ATG,observed_exact,nextclade_aligned_nt,3\n"
        "HA,2,2020,01,K,AAA,observed_exact,unmapped,2\n"
    )
    out = tmp_path / "codon_validation_report.csv"
    monkeypatch.setattr(
        "sys.argv",
        [
            "validate_codon_usage.py",
            "--count-root",
            str(count_root),
            "--out",
            str(out),
        ],
    )

    validate_main()

    source_summary = out.with_name("codon_validation_report_source_summary.csv").read_text()
    assert "HA,valid_observed_codon,nextclade_aligned_nt,3" in source_summary
    assert "HA,valid_observed_codon,unmapped,2" in source_summary
