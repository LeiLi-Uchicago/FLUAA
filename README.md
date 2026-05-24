# FLUAA

FLUAA is a Nextflow and Python pipeline for processing influenza nucleotide FASTA files and metadata by subtype. It prepares metadata, splits and deduplicates sequences by gene, runs Nextclade v3, annotates HA/NA clades back onto metadata, and generates amino acid plus codon count tables for downstream analysis.

Supported profiles:

| Profile | Subtype |
|---|---|
| `H1N1` | Influenza A/H1N1, split into pdm09 and seasonal references by metadata lineage |
| `H3N2` | Influenza A/H3N2 |
| `B_VIC` | Influenza B Victoria |
| `B_YAM` | Influenza B Yamagata |

## Input Layout (GISAID download format)

Each subtype input directory should contain nucleotide FASTA files and metadata Excel files.

Expected files:

```text
H1N1/
  *-NT.fasta
  *.xls or *.xlsx
```

FASTA headers should have this format:

```text
>Isolate_Name|Isolate_Id|gene
```

The metadata file must contain a unique isolate identifier column. By default the pipeline uses:

| Parameter | Default |
|---|---|
| `--id_column` | `Isolate_Id` |
| `--date_column` | `Collection_Date` |

Dates may be complete or partial. The pipeline derives `Year` and `Month` from the metadata date column.

## Requirements

The pipeline expects:

- Nextflow
- Python 3.11 or newer
- `nextclade3` on `PATH`, or available through the included conda environment
- Python packages listed in `requirements.txt`

Python dependencies:

```text
pandas
xlrd
openpyxl
biopython
pyarrow
pytest
```

A conda environment is provided at:

```text
envs/flu_pipeline.yml
```

The default command names are configured in `nextflow.config`:

| Parameter | Default |
|---|---|
| `--python_cmd` | `python` |
| `--nextclade_cmd` | `nextclade3` |
| `--python_env` | `envs/flu_pipeline.yml` |
| `--nextclade_env` | `envs/flu_pipeline.yml` |

If your environment uses `python3` instead of `python`, run with:

```bash
nextflow run main.nf -profile H1N1 --python_cmd python3
```

## Quick Start

Run one subtype with the default local input/output paths:

```bash
nextflow run main.nf -profile H1N1
```

Run with conda environment creation:

```bash
nextflow run main.nf -profile H1N1 -with-conda
```

Resume a previous run after fixing an error:

```bash
nextflow run main.nf -profile H1N1 -resume
```

Run a small smoke test:

```bash
nextflow run main.nf -profile H1N1 \
  --max_records_test 100
```

## Run Examples

These examples mirror `run.sh`.

Write results inside this repository:

```bash
nextflow run main.nf -profile H1N1 \
  --input_dir Data/H1N1 \
  --outdir results/H1N1
```

```bash
nextflow run main.nf -profile H3N2 \
  --input_dir Data/H3N2 \
  --outdir results/H3N2
```

```bash
nextflow run main.nf -profile B_VIC \
  --input_dir Data/B_VIC \
  --outdir results/B_VIC
```

```bash
nextflow run main.nf -profile B_YAM \
  --input_dir Data/B_YAM \
  --outdir results/B_YAM
```

## Main Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input_dir` | Profile-specific subtype folder | Input folder containing FASTA and metadata files |
| `--outdir` | `results/<subtype>` | Output directory |
| `--id_column` | `Isolate_Id` | Metadata column used as the unique isolate identifier |
| `--date_column` | `Collection_Date` | Metadata date column used to derive `Year` and `Month` |
| `--dedupe_mode` | `best_sequence` | Duplicate sequence handling mode |
| `--insertion_min_support` | `100` | Minimum strain support for insertion events included in supported insertion output |
| `--max_records_test` | `0` | Limit input records for smoke testing; `0` means no limit |
| `--python_cmd` | `python` | Python executable |
| `--nextclade_cmd` | `nextclade3` | Nextclade executable |

## Pipeline Steps

1. `PREPARE_INPUTS`

   Reads subtype metadata and FASTA files, normalizes metadata, writes merged metadata, splits nucleotide sequences by gene, and deduplicates records.

2. `BUILD_NEXTCLADE_MANIFEST`

   Builds the list of gene FASTA files and Nextclade datasets to run. For H1N1, records are split into `pdm09` and `seasonal` groups using the metadata `Lineage` column and year.

3. `RUN_NEXTCLADE`

   Runs Nextclade v3 for each gene/reference group and produces TSV, aligned nucleotide FASTA, translated protein FASTA, and an internal NDJSON file.

4. `PUBLISH_NEXTCLADE_RESULTS`

   Publishes compact Nextclade outputs to the final output folder. The large Nextclade NDJSON file is kept in Nextflow `work/` for insertion parsing but is not copied to `results/`.

5. `ANNOTATE_METADATA`

   Adds HA and NA clade annotations to the merged metadata and writes `metadata_merged_annotated.csv`.

6. `GENERATE_COUNTS`

   Generates amino acid and codon count tables. Codons are inferred only from Nextclade aligned nucleotide FASTA plus Nextclade translated amino acid FASTA. Raw input NT sequences are not used as a codon fallback.

7. `VALIDATE_CODONS`

   Validates amino acid/codon consistency and reports valid codons, mismatches, ambiguous codons, unrecovered codons, deletions, and source categories.

## Outputs

For each subtype output directory:

```text
results/H1N1/
  metadata_merged_annotated.csv
  codon_validation_report.csv
  codon_validation_report_summary.csv
  codon_validation_report_source_summary.csv
  prepared/
  nextclade/
  count/
```

Prepared input outputs:

```text
prepared/metadata/merged_metadata.csv
prepared/fasta_by_gene/<gene>.fasta
prepared/reports/prepare_summary.csv
prepared/reports/fasta_duplicate_conflicts.csv
prepared/reports/malformed_fasta_headers.csv
```

Published Nextclade outputs:

```text
nextclade/<group>_<gene>/nextclade.tsv
nextclade/<group>_<gene>/aligned.fasta
nextclade/<group>_<gene>/translations/<protein>.fasta
```

Count outputs:

```text
count/<protein>/aa_usage_by_Year_Month.csv
count/<protein>/aa_usage_by_<clade_column>.csv
count/insertions/insertion_events.csv
count/insertions/insertion_summary.csv
count/insertions/supported_insertions.csv
```

Count tables include:

| Column | Description |
|---|---|
| `Protein` | Protein or CDS name |
| `Position` | Amino acid position |
| `Year` | Year derived from metadata |
| `Month` | Month derived from metadata |
| `AminoAcid` | Amino acid state |
| `Codon` | Codon assigned to the amino acid state |
| `CodonStatus` | Codon classification, such as `observed_exact`, `codon_aa_mismatch`, `codon_ambiguous`, `deletion`, or `codon_unavailable` |
| `CodonSource` | Source used for codon recovery, usually `nextclade_aligned_nt`, `deletion`, or `unmapped` |
| `Count` | Aggregated observation count |

Validation outputs:

| File | Description |
|---|---|
| `codon_validation_report.csv` | Detailed row-level validation report |
| `codon_validation_report_summary.csv` | Counts by protein and validation issue |
| `codon_validation_report_source_summary.csv` | Counts by protein, validation issue, and codon source |

## Clade Groupings

The count tables are grouped by the clade columns configured per profile.

| Profile | Clade grouping columns |
|---|---|
| `H1N1` | `HA_clade`, `HA_short_clade`, `HA_legacy_clade`, `NA_clade` |
| `H3N2` | `HA_clade`, `HA_short_clade`, `HA_legacy_clade`, `NA_clade` |
| `B_VIC` | `HA_clade`, `HA_legacy_clade`, `NA_clade` |
| `B_YAM` | `HA_clade` |

## Nextclade References

The configured Nextclade references are summarized in `nextclade_references.md`.

| Subtype | Lineage / Group | HA | NA | MP | NP | NS | PA | PB1 | PB2 |
|---|---|---|---|---|---|---|---|---|---|
| H1N1 | pdm09 | `flu_h1n1pdm_ha` | `flu_h1n1pdm_na` | `nextstrain/flu/h1n1pdm/mp` | `nextstrain/flu/h1n1pdm/np` | `nextstrain/flu/h1n1pdm/ns` | `nextstrain/flu/h1n1pdm/pa` | `nextstrain/flu/h1n1pdm/pb1` | `nextstrain/flu/h1n1pdm/pb2` |
| H1N1 | seasonal | `flu_h1n1_ha` | `flu_h1n1_na` | `flu_h1n1_mp` | `flu_h1n1_np` | `flu_h1n1_ns` | `flu_h1n1_pa` | `flu_h1n1_pb1` | `flu_h1n1_pb2` |
| H3N2 | all | `nextstrain/flu/h3n2/ha/EPI1857216` | `nextstrain/flu/h3n2/na/EPI1857215` | `nextstrain/flu/h3n2/mp` | `nextstrain/flu/h3n2/np` | `nextstrain/flu/h3n2/ns` | `nextstrain/flu/h3n2/pa` | `nextstrain/flu/h3n2/pb1` | `nextstrain/flu/h3n2/pb2` |
| B_VIC | all | `nextstrain/flu/vic/ha/KX058884` | `nextstrain/flu/vic/na/CY073894` | `nextstrain/flu/vic/mp` | `nextstrain/flu/vic/np` | `nextstrain/flu/vic/ns` | `nextstrain/flu/vic/pa` | `nextstrain/flu/vic/pb1` | `nextstrain/flu/vic/pb2` |
| B_YAM | all | `nextstrain/flu/yam/ha/JN993010` | `nextstrain/flu/b/na/CY073894` | `nextstrain/flu/b/mp` | `nextstrain/flu/b/np` | `nextstrain/flu/b/ns` | `nextstrain/flu/b/pa` | `nextstrain/flu/b/pb1` | `nextstrain/flu/b/pb2` |

## H1N1 Lineage Split

The H1N1 profile uses `params.h1n1_split_lineage = true`.

The manifest builder assigns records as follows:

| Metadata condition | Nextclade group |
|---|---|
| `Lineage` is `pdm09` or a pdm09 alias | `pdm09` |
| `Lineage` is `seasonal` | `seasonal` |
| `Lineage` is empty and `Year < 2010` | `seasonal` |
| `Lineage` is empty and `Year >= 2010` or missing | `pdm09` |

