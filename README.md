# FLUAA

FLUAA is a Nextflow and Python pipeline for processing influenza nucleotide FASTA files and metadata by subtype. It prepares metadata, splits and deduplicates sequences by gene, runs Nextclade v3, annotates HA/NA clades back onto metadata, and generates amino acid plus codon count tables for downstream analysis.

Supported profiles:

| Profile | Subtype |
|---|---|
| `H1N1` | Influenza A/H1N1. Each isolate is classified from its HA sequence as **pdm09** or **seasonal**, mapped against its own reference set, and reported separately |
| `H3N2` | Influenza A/H3N2 |
| `B_VIC` | Influenza B Victoria |
| `B_YAM` | Influenza B Yamagata |
| `H5NX` | Avian influenza A/H5Nx, with HA clading and NA subtype-aware counts |

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

Dates may be complete or partial. The pipeline derives `Year` and `Month` from the metadata date column. Year-only and messy partial dates (e.g. `2009`, `2009.0`, `2009-XX-XX`, `2009 (Month unknown)`) are kept — the year is retained and the month is reported as `Unknown` rather than dropping the record. When the collection date is missing but the strain name ends in a plausible year (e.g. `A/USSR/90/1977`), that year is used; a trailing isolate number that is not a plausible year (e.g. `A/Penarth/6684`) is not mistaken for one.

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

Run the bundled H5Nx example:

```bash
nextflow run main.nf -profile H5NX \
  --input_dir h5_example \
  --outdir results/H5NX
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

```bash
nextflow run main.nf -profile H5NX \
  --input_dir Data/H5NX \
  --outdir results/H5NX
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
| `--h5_ha_dataset` | `community/moncla-lab/iav-h5/ha/all-clades` | H5 HA Nextclade dataset used by the `H5NX` profile |
| `--h1n1_split_lineage` | `true` (H1N1 profile) | Classify H1N1 isolates as pdm09/seasonal and process/report them separately. When `false`, all H1N1 is mapped against the pdm09 references with no split |
| `--h1n1_max_divergence_frac` | `0.15` | Max per-site HA nucleotide divergence for a lineage call to be trusted; isolates above it (matching neither reference) are filtered out |
| `--h1n1_pdm_min_year` | `2009` | pdm09-classified isolates collected before this year are filtered out (seasonal H1N1 stopped circulating after the 2009 pandemic) |

## Pipeline Steps

1. `PREPARE_INPUTS`

   Reads subtype metadata and FASTA files, normalizes metadata, writes merged metadata, splits nucleotide sequences by gene, and deduplicates records.

2. `CLASSIFY_H1N1_LINEAGE` *(H1N1 only, when `h1n1_split_lineage = true`)*

   Aligns every HA sequence to **both** the pdm09 reference and the custom seasonal reference with Nextclade, then assigns each isolate to whichever reference it matches with the lowest nucleotide divergence. A divergence gate and a pre-pandemic year filter drop isolates that match neither reference (see [H1N1 seasonal vs pdm09](#h1n1-seasonal-vs-pdm09)). Writes `h1n1_lineage.csv`, a detailed `h1n1_classification/` folder, and `lineage_filter_out.csv`.

3. `BUILD_NEXTCLADE_MANIFEST`

   Builds the list of gene FASTA files and Nextclade datasets to run. For H1N1 (split mode), records are split into `pdm09` and `seasonal` groups by the classifier's assignments and each group is routed to its matching reference set for all eight segments. For H5NX, NA records are split by parsed `NA_subtype`; only configured N1/N2 groups are sent to Nextclade and the rest are handled by the H5 NA fallback counter.

4. `RUN_NEXTCLADE`

   Runs Nextclade v3 for each gene/reference group and produces TSV, aligned nucleotide FASTA, translated protein FASTA, and an internal NDJSON file. Named community datasets are used with `--dataset-name`; local custom dataset directories (the seasonal H1N1 references under `refs/`) are used with `--input-ref` / `--input-annotation`.

5. `PUBLISH_NEXTCLADE_RESULTS`

   Publishes compact Nextclade outputs to the final output folder. The large Nextclade NDJSON file is kept in Nextflow `work/` for insertion parsing but is not copied to `results/`.

6. `ANNOTATE_METADATA`

   Adds HA and NA clade annotations to the merged metadata and writes `metadata_merged_annotated.csv`.

7. `GENERATE_COUNTS`

   Generates amino acid and codon count tables. Codons are inferred from Nextclade aligned nucleotide FASTA plus Nextclade translated amino acid FASTA. In H1N1 split mode, count tables are emitted per lineage (protein labels are suffixed `_pdm09` / `_seasonal`), the seasonal HA stalk deletion is normalized to a single position (see below), and seasonal tables are grouped by year-month only. For H5NX NA records that have no Nextclade translation, the pipeline infers an NA ORF from the raw nucleotide sequence and marks those codons with `CodonSource=fallback_orf`.

8. `ORGANIZE_BY_LINEAGE` *(H1N1 only, when `h1n1_split_lineage = true`)*

   Splits the annotated metadata and the per-lineage count tables into two self-contained folders, `H1N1pdm09/` and `H1N1seasonal/`, each with its own `metadata_merged_annotated.csv` and `count/` directory. Isolates classified `filtered_out` are excluded from both.

9. `VALIDATE_CODONS`

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

### H1N1 split-lineage outputs

When `h1n1_split_lineage = true` (the default for the H1N1 profile), the run additionally produces the lineage classification and two self-contained per-lineage result folders:

```text
results/H1N1/
  h1n1_lineage.csv                       # Isolate_Id -> pdm09 | seasonal | filtered_out (+ divergence, confidence)
  lineage_filter_out.csv                 # isolates dropped, with full metadata and the reason
  h1n1_classification/
    divergence_details.csv               # per-isolate divergence to both references
    classification_summary.txt           # counts, confidence distribution, thresholds used
  H1N1pdm09/
    metadata_merged_annotated.csv
    count/<protein>/...
  H1N1seasonal/
    metadata_merged_annotated.csv
    count/<protein>/...
```

The unsplit `count/` directory still holds the combined tables (with `_pdm09` / `_seasonal` protein suffixes); `H1N1pdm09/` and `H1N1seasonal/` are the tidied, per-lineage views with the suffix removed.

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
| `H1N1` (pdm09) | `HA_clade`, `HA_short_clade`, `HA_legacy_clade`, `NA_clade` |
| `H1N1` (seasonal) | none — year-month only (no maintained clade scheme for extinct seasonal H1N1) |
| `H3N2` | `HA_clade`, `HA_short_clade`, `HA_legacy_clade`, `NA_clade` |
| `B_VIC` | `HA_clade`, `HA_legacy_clade`, `NA_clade` |
| `B_YAM` | `HA_clade` |
| `H5NX` | `HA_clade`, `NA_subtype`, `Pathogenicity` |

## Nextclade References

| Subtype | Lineage / Group | HA | NA | MP | NP | NS | PA | PB1 | PB2 |
|---|---|---|---|---|---|---|---|---|---|
| H1N1 | pdm09 | `flu_h1n1pdm_ha` | `flu_h1n1pdm_na` | `nextstrain/flu/h1n1pdm/mp` | `nextstrain/flu/h1n1pdm/np` | `nextstrain/flu/h1n1pdm/ns` | `nextstrain/flu/h1n1pdm/pa` | `nextstrain/flu/h1n1pdm/pb1` | `nextstrain/flu/h1n1pdm/pb2` |
| H1N1 | seasonal | `refs/ussr77_h1n1_ha` | `refs/ussr77_h1n1_na` | `refs/ussr77_h1n1_mp` | `refs/ussr77_h1n1_np` | `refs/ussr77_h1n1_ns` | `refs/ussr77_h1n1_pa` | `refs/ussr77_h1n1_pb1` | `refs/ussr77_h1n1_pb2` |
| H3N2 | all | `nextstrain/flu/h3n2/ha/EPI1857216` | `nextstrain/flu/h3n2/na/EPI1857215` | `nextstrain/flu/h3n2/mp` | `nextstrain/flu/h3n2/np` | `nextstrain/flu/h3n2/ns` | `nextstrain/flu/h3n2/pa` | `nextstrain/flu/h3n2/pb1` | `nextstrain/flu/h3n2/pb2` |
| B_VIC | all | `nextstrain/flu/vic/ha/KX058884` | `nextstrain/flu/vic/na/CY073894` | `nextstrain/flu/vic/mp` | `nextstrain/flu/vic/np` | `nextstrain/flu/vic/ns` | `nextstrain/flu/vic/pa` | `nextstrain/flu/vic/pb1` | `nextstrain/flu/vic/pb2` |
| B_YAM | all | `nextstrain/flu/yam/ha/JN993010` | `nextstrain/flu/b/na/CY073894` | `nextstrain/flu/b/mp` | `nextstrain/flu/b/np` | `nextstrain/flu/b/ns` | `nextstrain/flu/b/pa` | `nextstrain/flu/b/pb1` | `nextstrain/flu/b/pb2` |
| H5NX | HA/all + internals | `community/moncla-lab/iav-h5/ha/all-clades` | N1: `nextstrain/flu/h1n1pdm/na/MW626056`; N2: `nextstrain/flu/h2n2/na`; N3-N9: fallback ORF | `nextstrain/flu/h2n2/mp` | `nextstrain/flu/h2n2/np` | `nextstrain/flu/h2n2/ns` | `nextstrain/flu/h2n2/pa` | `nextstrain/flu/h2n2/pb1` | `nextstrain/flu/h2n2/pb2` |

## H5Nx Support

The `H5NX` profile uses the official H5 HA Nextclade all-clades dataset by default. You can override it for focused HA clade runs:

```bash
nextflow run main.nf -profile H5NX \
  --h5_ha_dataset community/moncla-lab/iav-h5/ha/2.3.4.4
```

The current H5 Nextclade datasets are HA-only, so FLUAA uses H2N2 influenza A datasets for internal-gene alignment/translation and subtype-aware NA outputs. NA count tables are named `NA_N1` through `NA_N9`. N1/N2 can use Nextclade translations when configured; other NA subtypes use fallback ORF translation with `CodonSource=fallback_orf`. NA records whose subtype cannot be parsed are skipped from NA protein counts.

## H1N1 seasonal vs pdm09

Pre-2009 seasonal H1N1 and pandemic H1N1pdm09 are distinct lineages with
different HA/NA reference sequences and **incompatible position numbering**.
Forcing both onto one reference produces alignment artifacts (spurious deletions
and unalignable codons), so the H1N1 profile keeps them apart end to end. This is
controlled by `params.h1n1_split_lineage = true` (the default for the profile).

### Classification (sequence-based)

Each isolate is classified from its **HA sequence**, not from metadata (GISAID
lineage/date labels are sometimes wrong). `CLASSIFY_H1N1_LINEAGE` aligns HA to
both references and assigns the isolate to whichever it matches with lower
nucleotide divergence (`totalSubstitutions + totalDeletions + totalInsertions`):

- A **failed alignment** to one reference scores as maximally divergent
  (`no_alignment`), never as a perfect match — so a seasonal strain that cannot
  align to the pdm09 reference is not mistaken for a perfect pdm09 match.
- A **divergence gate** (`--h1n1_max_divergence_frac`, default `0.15` per site,
  normalized by aligned length) rejects "wins" that are only close relative to
  the ~22% pdm09-vs-seasonal distance. Isolates matching neither reference are
  written to `lineage_filter_out.csv` as `filtered_out`.
- A **pre-pandemic filter** (`--h1n1_pdm_min_year`, default `2009`) drops any
  isolate classified pdm09 but collected before that year, since seasonal H1N1
  stopped circulating in humans after 2009.

Outputs of this step:

| File | Description |
|---|---|
| `h1n1_lineage.csv` | Per-isolate `H1_lineage` (`pdm09` / `seasonal` / `filtered_out`), `confidence`, divergence scores, and `gap_size` |
| `h1n1_classification/divergence_details.csv` | Full per-isolate divergence to both references (absolute and per-site) |
| `h1n1_classification/classification_summary.txt` | Totals, confidence distribution, and thresholds used |
| `lineage_filter_out.csv` | Filtered-out isolates with full metadata and the drop reason |

### Seasonal references

Seasonal genes are aligned to custom Nextclade datasets built from
**A/USSR/90/1977 (EPI_ISL_66104)** for all eight segments, stored as local dataset
directories under `refs/ussr77_h1n1_<gene>/` (each with `reference.fasta`,
`genome_annotation.gff3`, and `pathogen.json`). These are alignment/translation
only — no tree or clade assignment — which is why seasonal count tables are
grouped by year-month rather than clade.

### Seasonal HA gap standardization

The seasonal HA stalk deletion can be placed at slightly different positions by
the aligner from strain to strain. `GENERATE_COUNTS` normalizes any gap within a
window of the canonical position (default **147**, counted from the start Met) to
that single position, so the same biological deletion is counted consistently
across all seasonal strains.

### Disabling the split

Run with `--h1n1_split_lineage false` to map all H1N1 against the pdm09
references with no classification or splitting (the classifier and organize steps
are skipped). This reproduces the older single-reference behavior.
