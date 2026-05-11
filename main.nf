#!/usr/bin/env nextflow

import groovy.json.JsonOutput

nextflow.enable.dsl = 2

def requiredParams = ['subtype', 'input_dir', 'datasets']
requiredParams.each { key ->
  if (!params.containsKey(key) || params[key] == null) {
    error "Missing params.${key}. Run with one of: -profile H1N1, H3N2, B_VIC, B_YAM"
  }
}

process PREPARE_INPUTS {
  tag "${subtype}"
  publishDir "${params.outdir}/prepared", mode: 'copy'
  conda params.python_env

  input:
  val subtype
  val input_dir

  output:
  path "prepared/metadata/merged_metadata.csv", emit: metadata
  path "prepared/fasta_by_gene/*.fasta", emit: gene_fastas
  path "prepared/reports", emit: reports

  script:
  def maxOpt = params.max_records_test as Integer
  """
  export PYTHONPATH='${projectDir}':\${PYTHONPATH:-}
  ${params.python_cmd} ${projectDir}/scripts/prepare_inputs.py \\
    --input-dir '${input_dir}' \\
    --outdir prepared \\
    --subtype '${subtype}' \\
    --id-column '${params.id_column}' \\
    --date-column '${params.date_column}' \\
    --dedupe-mode '${params.dedupe_mode}' \\
    --max-records-test ${maxOpt}
  """

  stub:
  """
  mkdir -p prepared/metadata prepared/fasta_by_gene prepared/reports
  printf 'Isolate_Id,Isolate_Name,strain,Collection_Date,Year,Month,Subtype,Lineage\\nEPI_STUB,stub|EPI_STUB|HA,stub,2020-01-01,2020,01,${subtype},stub\\n' > prepared/metadata/merged_metadata.csv
  printf '>stub|EPI_STUB|HA\\nATGAAA\\n' > prepared/fasta_by_gene/HA.fasta
  printf '>stub|EPI_STUB|NA\\nATGAAA\\n' > prepared/fasta_by_gene/NA.fasta
  printf 'metric,value\\nstub,1\\n' > prepared/reports/prepare_summary.csv
  : > prepared/reports/fasta_duplicate_conflicts.csv
  : > prepared/reports/malformed_fasta_headers.csv
  """
}

process RUN_NEXTCLADE {
  tag "${group}:${gene}"
  conda params.nextclade_env

  input:
  tuple val(group), val(gene), path(fasta), val(dataset)

  output:
  path "${group}_${gene}", emit: result_dirs

  script:
  """
  mkdir -p '${group}_${gene}/translations'
  ${params.nextclade_cmd} run \\
    --dataset-name '${dataset}' \\
    --output-tsv '${group}_${gene}/nextclade.tsv' \\
    --output-ndjson '${group}_${gene}/nextclade.ndjson' \\
    --output-fasta '${group}_${gene}/aligned.fasta' \\
    --output-translations '${group}_${gene}/translations/{cds}.fasta' \\
    '${fasta}'
  """

  stub:
  """
  mkdir -p '${group}_${gene}/translations'
  printf 'seqName\\tclade\\tproposedSubclade\\tsubclade\\tshort_clade\\tlegacy_clade\\n' > '${group}_${gene}/nextclade.tsv'
  printf 'stub|EPI_STUB|${gene}\\tstub-clade\\tstub-proposed\\tstub-sub\\tstub-short\\tstub-legacy\\n' >> '${group}_${gene}/nextclade.tsv'
  printf '{"seqName":"stub|EPI_STUB|${gene}"}\\n' > '${group}_${gene}/nextclade.ndjson'
  printf '>stub|EPI_STUB|${gene}\\nATGAAA\\n' > '${group}_${gene}/aligned.fasta'
  printf '>stub|EPI_STUB|${gene}\\nMK\\n' > '${group}_${gene}/translations/${gene}.fasta'
  """
}

process PUBLISH_NEXTCLADE_RESULTS {
  tag "${params.subtype}"
  publishDir "${params.outdir}/nextclade", mode: 'copy', saveAs: { filename ->
    filename.startsWith('nextclade_published/') ? filename.substring('nextclade_published/'.length()) : filename
  }

  input:
  path nextclade_dirs

  output:
  path "nextclade_published/*", emit: published_dirs

  script:
  def dirs = nextclade_dirs.collect { "'${it}'" }.join(' ')
  """
  mkdir -p nextclade_published
  for dir in ${dirs}; do
    out="nextclade_published/\$(basename "\$dir")"
    mkdir -p "\$out/translations"
    cp "\$dir/nextclade.tsv" "\$out/"
    cp "\$dir/aligned.fasta" "\$out/"
    cp "\$dir"/translations/*.fasta "\$out/translations/" 2>/dev/null || true
  done
  """

  stub:
  """
  mkdir -p nextclade_published/stub/translations
  printf 'seqName\\tclade\\n' > nextclade_published/stub/nextclade.tsv
  printf '>stub|EPI_STUB|HA\\nATGAAA\\n' > nextclade_published/stub/aligned.fasta
  printf '>stub|EPI_STUB|HA\\nMK\\n' > nextclade_published/stub/translations/HA.fasta
  """
}

process BUILD_NEXTCLADE_MANIFEST {
  tag "${params.subtype}"
  conda params.python_env

  input:
  path metadata
  path gene_fastas

  output:
  path "nextclade_manifest.csv", emit: manifest
  path "split_fastas", emit: split_fastas

  script:
  def fastas = gene_fastas.collect { "'${it}'" }.join(' ')
  def datasetsJson = JsonOutput.toJson(params.datasets ?: [:])
  def seasonalJson = JsonOutput.toJson(params.h1n1_seasonal_datasets ?: [:])
  def splitFlag = params.h1n1_split_lineage ? "--h1n1-split-lineage" : ""
  """
  export PYTHONPATH='${projectDir}':\${PYTHONPATH:-}
  ${params.python_cmd} ${projectDir}/scripts/build_nextclade_manifest.py \\
    --metadata '${metadata}' \\
    --gene-fastas ${fastas} \\
    --outdir . \\
    --subtype '${params.subtype}' \\
    --id-column '${params.id_column}' \\
    --datasets-json '${datasetsJson}' \\
    --h1n1-seasonal-datasets-json '${seasonalJson}' \\
    ${splitFlag}
  """

  stub:
  """
  mkdir -p split_fastas/all
  printf '>stub|EPI_STUB|HA\\nATGAAA\\n' > split_fastas/all/HA.fasta
  printf '>stub|EPI_STUB|NA\\nATGAAA\\n' > split_fastas/all/NA.fasta
  printf 'group,gene,fasta,dataset\\n' > nextclade_manifest.csv
  printf 'all,HA,split_fastas/all/HA.fasta,${(params.datasets ?: [HA:"stub"]).HA ?: "stub"}\\n' >> nextclade_manifest.csv
  printf 'all,NA,split_fastas/all/NA.fasta,${(params.datasets ?: [NA:"stub"]).NA ?: "stub"}\\n' >> nextclade_manifest.csv
  """
}

process ANNOTATE_METADATA {
  tag "${params.subtype}"
  publishDir "${params.outdir}", mode: 'copy'
  conda params.python_env

  input:
  path metadata
  path nextclade_dirs

  output:
  path "metadata_merged_annotated.csv", emit: annotated_metadata

  script:
  def dirs = nextclade_dirs.collect { "'${it}'" }.join(' ')
  """
  export PYTHONPATH='${projectDir}':\${PYTHONPATH:-}
  ${params.python_cmd} ${projectDir}/scripts/annotate_metadata.py \\
    --metadata '${metadata}' \\
    --subtype '${params.subtype}' \\
    --out metadata_merged_annotated.csv \\
    --nextclade-dirs ${dirs}
  """

  stub:
  """
  printf 'Isolate_Id,Isolate_Name,strain,Collection_Date,Year,Month,Subtype,Lineage,HA_clade,HA_proposedSubclade,HA_subclade,HA_short_clade,HA_legacy_clade,NA_clade\\n' > metadata_merged_annotated.csv
  printf 'EPI_STUB,stub|EPI_STUB|HA,stub,2020-01-01,2020,01,${params.subtype},stub,stub-clade,stub-proposed,stub-sub,stub-short,stub-legacy,stub-clade\\n' >> metadata_merged_annotated.csv
  """
}

process GENERATE_COUNTS {
  tag "${params.subtype}"
  publishDir "${params.outdir}", mode: 'copy'
  conda params.python_env

  input:
  path annotated_metadata
  path nextclade_dirs
  path gene_fastas

  output:
  path "count", emit: count_outputs

  script:
  def dirs = nextclade_dirs.collect { "'${it}'" }.join(' ')
  def fastas = gene_fastas.collect { "'${it}'" }.join(' ')
  def clades = (params.clade_columns ?: []).join(',')
  """
  export PYTHONPATH='${projectDir}':\${PYTHONPATH:-}
  echo 'count_code_version: nextclade_aligned_nt_only_v7' >&2
  ${params.python_cmd} ${projectDir}/scripts/generate_counts.py \\
    --metadata '${annotated_metadata}' \\
    --nextclade-dirs ${dirs} \\
    --gene-fastas ${fastas} \\
    --outdir count \\
    --insertion-min-support ${params.insertion_min_support} \\
    --clade-columns '${clades}'
  """

  stub:
  """
  mkdir -p count/stub count/insertions
  printf 'Protein,Position,Year,Month,AminoAcid,Codon,CodonStatus,CodonSource,Count\\n' > count/stub/aa_usage_by_Year_Month.csv
  printf 'seq_name,Isolate_Id,protein,position_label,aa_state,codon,Year,Month\\n' > count/insertions/insertion_events.csv
  printf 'protein,position_label,aa_state,codon,support\\n' > count/insertions/insertion_summary.csv
  printf 'protein,position_label,aa_state,codon,support\\n' > count/insertions/supported_insertions.csv
  """
}

process VALIDATE_CODONS {
  tag "${params.subtype}"
  publishDir "${params.outdir}", mode: 'copy'
  conda params.python_env

  input:
  path count_dir

  output:
  path "codon_validation_report.csv", emit: detailed_report
  path "codon_validation_report_summary.csv", emit: summary_report
  path "codon_validation_report_source_summary.csv", emit: source_summary_report

  script:
  """
  export PYTHONPATH='${projectDir}':\${PYTHONPATH:-}
  ${params.python_cmd} ${projectDir}/scripts/validate_codon_usage.py \\
    --count-root '${count_dir}' \\
    --out codon_validation_report.csv
  """

  stub:
  """
  printf 'Protein,issue,Position,Grouping,GroupingValue,Year,Month,AminoAcid,Codon,CodonStatus,CodonSource,TranslatedAA,Count,SourceTable\\n' > codon_validation_report.csv
  printf 'Protein,issue,Count\\n' > codon_validation_report_summary.csv
  printf 'Protein,issue,CodonSource,Count\\n' > codon_validation_report_source_summary.csv
  """
}

workflow {
  prepared = PREPARE_INPUTS(params.subtype, params.input_dir)

  prepared_gene_fastas = prepared.gene_fastas.collect()
  manifest = BUILD_NEXTCLADE_MANIFEST(prepared.metadata, prepared_gene_fastas)

  nextclade_inputs = manifest.manifest
    .splitCsv(header: true)
    .map { row -> tuple(row.group.toString(), row.gene.toString().toUpperCase(), file(row.fasta.toString()), row.dataset.toString()) }
  nextclade_results = RUN_NEXTCLADE(nextclade_inputs)
  nextclade_dirs = nextclade_results.result_dirs.collect()
  PUBLISH_NEXTCLADE_RESULTS(nextclade_dirs)

  annotated = ANNOTATE_METADATA(prepared.metadata, nextclade_dirs)
  counts = GENERATE_COUNTS(annotated.annotated_metadata, nextclade_dirs, prepared_gene_fastas)
  VALIDATE_CODONS(counts.count_outputs)
}
