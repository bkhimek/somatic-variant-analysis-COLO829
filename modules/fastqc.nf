// Module 1 — QC
// FastQC per sample (tumour and normal run through this independently, per plan §6),
// then MultiQC aggregates both into one report — this cross-sample aggregation is the
// one thing Phase 1 adds beyond a straight reuse of Project 4's QC step.

process FASTQC {
    // tag/publishDir written as closures (found via execution 2026-08-30, Nextflow 26.x):
    // a plain GString directive like tag "${sample_id}" is evaluated before task inputs are
    // bound, so sample_id isn't in scope yet -- "No such variable: sample_id". Older Nextflow
    // (24.10.5) silently made this work via implicit lazy-evaluation of GStrings referencing
    // input names; 26.x no longer does, so it must be an explicit closure. Same fix applied to
    // every other tag/publishDir in modules/alignment.nf that references an input variable.
    tag { sample_id }
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir { "${params.outdir}/qc/fastqc/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(reads1), path(reads2)

    output:
    tuple val(sample_id), path("*.html"), path("*.zip"), emit: reports

    script:
    """
    fastqc --threads ${task.cpus} ${reads1} ${reads2}
    """
}

process MULTIQC {
    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'
    publishDir "${params.outdir}/qc/multiqc", mode: 'copy'

    input:
    path('fastqc_reports/*')

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_report_data"), emit: data

    script:
    """
    multiqc fastqc_reports --filename multiqc_report.html
    """
    // NOTE (found via real execution, 2026-08-29): MultiQC names its data directory from the
    // --filename stem, not a fixed "multiqc_data" -- "--filename multiqc_report.html" produces
    // "multiqc_report_data", confirmed by an actual run's "Data : multiqc_report_data" log line.
    // The original "multiqc_data" guess (never executed, only manually reviewed) was wrong.
}
