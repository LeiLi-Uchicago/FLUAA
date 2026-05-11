# FLUAA Flu Nextflow Pipeline

This repository contains a Nextflow + Python pipeline for subtype-specific influenza processing.

## Run

If your current environment already has `nextclade3`, `xlrd`, `pandas`, and `biopython` available, run directly:

```bash
nextflow run main.nf -profile H1N1
```

Or create the conda environment automatically through Nextflow:

```bash
nextflow run main.nf -profile H1N1 -with-conda
```

Available profiles:

```bash
H1N1
H3N2
B_VIC
B_YAM
```

Useful parameters:

```bash
nextflow run main.nf -profile H3N2 -with-conda \
  --outdir results \
  --insertion_min_support 2 \
  --max_records_test 100
```

The Nextclade executable defaults to `nextclade3`; override with `--nextclade_cmd nextclade` if needed.
The Python executable defaults to `python`; override with `--python_cmd python3` or a full path if needed.

This config passes `--solver classic` to `conda env create`, which avoids the common broken `conda-libmamba-solver` / `libmambapy` error.

`--max_records_test 100` is intended for smoke tests and limits the number of FASTA records prepared.

## Outputs

For each profile, the pipeline writes:

- `prepared/metadata/merged_metadata.csv`
- `prepared/fasta_by_gene/*.fasta`
- `metadata_merged_annotated.csv`
- `count/<protein>/aa_usage_by_Year_Month.csv`
- `count/<protein>/aa_usage_by_<clade_column>.csv`
- `count/insertions/insertion_events.csv`
- `count/insertions/insertion_summary.csv`
- `count/insertions/supported_insertions.csv`
- `codon_validation_report.csv`
- `codon_validation_report_summary.csv`

Nextclade outputs are published under `nextclade/`.

## Notes

- H1N1 uses H1N1pdm Nextclade datasets for all records.
- B_YAM HA uses the Yamagata HA dataset, while B_YAM NA uses the all-influenza-B NA dataset.
- Nextclade v3 strips insertions from translated FASTA outputs, so insertion events are read separately from the NDJSON output and added to the main count tables only when they meet `--insertion_min_support`.
