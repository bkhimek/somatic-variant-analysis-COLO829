// Module 1 — QC
// FastQC per sample (tumour and normal run through this independently, per plan §6),
// then MultiQC aggregates both into one report — this cross-sample aggregation is the
// one thing Phase 1 adds beyond a straight reuse of Project 4's QC step.

process FASTQC {
    tag "${sample_id}"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/qc/fastqc/${sample_id}", mode: 'copy'

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
    path("multiqc_data"), emit: data

    script:
    """
    multiqc fastqc_reports --filename multiqc_report.html
    """
}
