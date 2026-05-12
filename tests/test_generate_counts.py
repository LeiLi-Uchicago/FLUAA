from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from scripts.generate_counts import (
    choose_codon_mapping,
    classify_aa,
    codons_for_aligned_peptide,
    main as counts_main,
    map_codons_from_nt_to_gapped_aa,
)


def test_classify_aa_states() -> None:
    assert classify_aa("A", "GCT") == ("A", "observed", "GCT")
    assert classify_aa("-", "AAA") == ("-", "deletion", "DEL")
    assert classify_aa("X", "NNN") == ("X", "ambiguous", "NA")
    assert classify_aa("?", "AAA") == ("?", "missing", "NA")


def test_aligned_peptide_deletions_do_not_consume_codons() -> None:
    assert codons_for_aligned_peptide("M-K", ["ATG", "AAA"]) == ["ATG", "DEL", "AAA"]


def test_direct_nt_to_gapped_aa_mapping_with_utr_and_ambiguous_codon() -> None:
    codons = map_codons_from_nt_to_gapped_aa("GGGCCCATGAAAACANNNTTTCCC", "---MKT-XF--")
    assert codons == ["DEL", "DEL", "DEL", "ATG", "AAA", "ACA", "DEL", "NNN", "TTT", "DEL", "DEL"]


def test_direct_nt_to_aa_mapping_can_skip_translated_intronic_region() -> None:
    codons = map_codons_from_nt_to_gapped_aa("ATGAAA" + "GCT" * 20 + "CCC", "MKP")
    assert codons == ["ATG", "AAA", "CCC"]


def test_codon_mapping_uses_only_nextclade_aligned_nt() -> None:
    codons = choose_codon_mapping("MK", aligned_nt="ATGAAA")
    assert codons == [("ATG", "nextclade_aligned_nt"), ("AAA", "nextclade_aligned_nt")]


def test_codon_mapping_without_aligned_nt_is_unmapped() -> None:
    codons = choose_codon_mapping("M-", aligned_nt=None)
    assert codons == [("NA", "unmapped"), ("DEL", "unmapped")]


def test_generate_counts_outputs_usage_and_insertions(tmp_path: Path, monkeypatch) -> None:
    metadata = tmp_path / "metadata.csv"
    pd.DataFrame(
        [
            {"Isolate_Id": "EPI_ISL_1", "Year": "2020", "Month": "07", "HA_clade": "6B"},
            {"Isolate_Id": "EPI_ISL_2", "Year": "2020", "Month": "08", "HA_clade": "6B"},
        ]
    ).to_csv(metadata, index=False)
    ha_dir = tmp_path / "HA"
    translations = ha_dir / "translations"
    translations.mkdir(parents=True)
    (ha_dir / "aligned.fasta").write_text(
        ">A/Test/1/2020|EPI_ISL_1|HA\nATGAAA\n>A/Test/2/2020|EPI_ISL_2|HA\nATG---\n"
    )
    (translations / "HA.fasta").write_text(
        ">A/Test/1/2020|EPI_ISL_1|HA\nMK\n>A/Test/2/2020|EPI_ISL_2|HA\nM-\n"
    )
    (ha_dir / "nextclade.ndjson").write_text(
        json.dumps(
            {
                "seqName": "A/Test/1/2020|EPI_ISL_1|HA",
                "annotation": annotation_for_ranges({"HA": (0, 6)}),
                "insertions": [{"cdsName": "HA", "pos": 2, "aa": "K", "codon": "AAA"}],
            }
        )
        + "\n"
        + json.dumps(
            {
                "seqName": "A/Test/2/2020|EPI_ISL_2|HA",
                "annotation": annotation_for_ranges({"HA": (0, 6)}),
            }
        )
        + "\n"
    )
    outdir = tmp_path / "counts"
    monkeypatch.setattr(
        "sys.argv",
        [
            "generate_counts.py",
            "--metadata",
            str(metadata),
            "--nextclade-dirs",
            str(ha_dir),
            "--outdir",
            str(outdir),
            "--insertion-min-support",
            "1",
            "--clade-columns",
            "HA_clade",
        ],
    )
    counts_main()
    year_month = outdir / "HA" / "aa_usage_by_Year_Month.csv"
    by_clade = outdir / "HA" / "aa_usage_by_HA_clade.csv"
    assert year_month.exists()
    assert by_clade.exists()
    assert year_month.read_text().splitlines()[0] == "Protein,Position,Year,Month,AminoAcid,Codon,CodonStatus,CodonSource,Count"
    assert by_clade.read_text().splitlines()[0] == "Protein,Position,HA_clade,Year,Month,AminoAcid,Codon,CodonStatus,CodonSource,Count"
    text = year_month.read_text()
    assert "HA,1,2020,07,M,ATG,observed_exact,nextclade_aligned_nt,1" in text
    assert "HA,2,2020,08,-,DEL,deletion,deletion,1" in text
    assert "insertion" not in text
    assert (outdir / "insertions" / "insertion_summary.csv").exists()
    assert (outdir / "insertions" / "supported_insertions.csv").exists()


def test_generate_counts_joins_split_ha_and_skips_pax(tmp_path: Path, monkeypatch) -> None:
    metadata = tmp_path / "metadata.csv"
    pd.DataFrame([{"Isolate_Id": "EPI_ISL_1", "Year": "2020", "Month": "", "HA_clade": ""}]).to_csv(metadata, index=False)
    ha_dir = tmp_path / "HA"
    ha_translations = ha_dir / "translations"
    ha_translations.mkdir(parents=True)
    (ha_dir / "aligned.fasta").write_text(">A/Test/1/2020|EPI_ISL_1|HA\nNNNNNNATGAAACCCGGG\n")
    (ha_translations / "SigPep.fasta").write_text(">A/Test/1/2020|EPI_ISL_1|HA\nM\n")
    (ha_translations / "HA1.fasta").write_text(">A/Test/1/2020|EPI_ISL_1|HA\nKP\n")
    (ha_translations / "HA2.fasta").write_text(">A/Test/1/2020|EPI_ISL_1|HA\nG\n")
    (ha_dir / "nextclade.ndjson").write_text(
        json.dumps(
            {
                "seqName": "A/Test/1/2020|EPI_ISL_1|HA",
                "annotation": annotation_for_ranges({"SigPep": (6, 9), "HA1": (9, 15), "HA2": (15, 18)}),
            }
        )
        + "\n"
    )

    pa_dir = tmp_path / "PA"
    pa_translations = pa_dir / "translations"
    pa_translations.mkdir(parents=True)
    (pa_dir / "aligned.fasta").write_text(">A/Test/1/2020|EPI_ISL_1|PA\nATGAAA\n")
    (pa_translations / "PA-X.fasta").write_text(">A/Test/1/2020|EPI_ISL_1|PA\nMK\n")
    (pa_translations / "PA.fasta").write_text(">A/Test/1/2020|EPI_ISL_1|PA\nMK\n")
    (pa_dir / "nextclade.ndjson").write_text(
        json.dumps(
            {
                "seqName": "A/Test/1/2020|EPI_ISL_1|PA",
                "annotation": annotation_for_ranges({"PA": (0, 6)}),
            }
        )
        + "\n"
    )

    outdir = tmp_path / "counts"
    monkeypatch.setattr(
        "sys.argv",
        [
            "generate_counts.py",
            "--metadata",
            str(metadata),
            "--nextclade-dirs",
            str(ha_dir),
            str(pa_dir),
            "--outdir",
            str(outdir),
            "--clade-columns",
            "HA_clade",
        ],
    )
    counts_main()

    text = (outdir / "HA" / "aa_usage_by_Year_Month.csv").read_text()
    assert "HA,1,2020,Unknown,M,ATG,observed_exact,nextclade_aligned_nt,1" in text
    assert "HA,2,2020,Unknown,K,AAA,observed_exact,nextclade_aligned_nt,1" in text
    assert "HA,3,2020,Unknown,P,CCC,observed_exact,nextclade_aligned_nt,1" in text
    assert "HA,4,2020,Unknown,G,GGG,observed_exact,nextclade_aligned_nt,1" in text
    assert (outdir / "PA" / "aa_usage_by_Year_Month.csv").exists()
    assert not (outdir / "PA-X").exists()


def test_generate_counts_merges_same_protein_from_multiple_nextclade_groups(tmp_path: Path, monkeypatch) -> None:
    metadata = tmp_path / "metadata.csv"
    pd.DataFrame(
        [
            {"Isolate_Id": "EPI_OLD", "Year": "2008", "Month": "01", "HA_clade": "seasonal"},
            {"Isolate_Id": "EPI_NEW", "Year": "2020", "Month": "01", "HA_clade": "pdm09"},
        ]
    ).to_csv(metadata, index=False)

    seasonal = tmp_path / "seasonal_HA"
    seasonal_translations = seasonal / "translations"
    seasonal_translations.mkdir(parents=True)
    (seasonal / "aligned.fasta").write_text(">old|EPI_OLD|HA\nATGAAA\n")
    (seasonal / "nextclade.ndjson").write_text('{"seqName":"old|EPI_OLD|HA"}\n')
    (seasonal_translations / "HA.fasta").write_text(">old|EPI_OLD|HA\nMK\n")

    pdm09 = tmp_path / "pdm09_HA"
    pdm09_translations = pdm09 / "translations"
    pdm09_translations.mkdir(parents=True)
    (pdm09 / "aligned.fasta").write_text(">new|EPI_NEW|HA\nATGCCC\n")
    (pdm09 / "nextclade.ndjson").write_text('{"seqName":"new|EPI_NEW|HA"}\n')
    (pdm09_translations / "HA.fasta").write_text(">new|EPI_NEW|HA\nMP\n")

    outdir = tmp_path / "counts"
    monkeypatch.setattr(
        "sys.argv",
        [
            "generate_counts.py",
            "--metadata",
            str(metadata),
            "--nextclade-dirs",
            str(seasonal),
            str(pdm09),
            "--outdir",
            str(outdir),
            "--clade-columns",
            "HA_clade",
        ],
    )
    counts_main()

    text = (outdir / "HA" / "aa_usage_by_Year_Month.csv").read_text()
    assert "HA,2,2008,01,K,AAA,observed_exact,nextclade_aligned_nt,1" in text
    assert "HA,2,2020,01,P,CCC,observed_exact,nextclade_aligned_nt,1" in text


def annotation_for_ranges(ranges: dict[str, tuple[int, int]]) -> dict:
    genes = []
    for name, (begin, end) in ranges.items():
        genes.append(
            {
                "name": name,
                "cdses": [
                    {
                        "segments": [
                            {
                                "range": {
                                    "begin": begin,
                                    "end": end,
                                },
                                "phase": 0,
                            }
                        ],
                        "proteins": [
                            {
                                "segments": [
                                    {
                                        "range": {
                                            "begin": begin,
                                            "end": end,
                                        }
                                    }
                                ]
                            }
                        ]
                    }
                ],
            }
        )
    return {"genes": genes}
