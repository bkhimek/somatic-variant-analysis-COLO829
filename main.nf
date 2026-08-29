#!/usr/bin/env nextflow
/*
 * COLO829 Tumour-Normal Somatic Genomics Pipeline
 * Entry point — Phase 1 wires up Modules 1-2 (QC, Alignment) via workflows/somatic.nf.
 * Phase 2+ will extend the same SOMATIC workflow rather than adding new entry points.
 *
 * Run:
 *   conda activate nextflow
 *   nextflow run main.nf -profile docker,dev  -params-file conf/dev_params.yaml   # fast smoke test
 *   nextflow run main.nf -profile docker,full -params-file conf/full_params.yaml  # real WGS run
 * (conf/*_params.yaml are not created yet -- see docs/PHASE1_NOTES.md; until then, pass
 * every required param on the command line as shown in the run instructions.)
 */

nextflow.enable.dsl = 2

include { SOMATIC } from './workflows/somatic.nf'

// ---- Required params: fail fast with a clear message rather than a cryptic NPE mid-DAG ----
def required = [
    'tumour_reads_1', 'tumour_reads_2',
    'normal_reads_1', 'normal_reads_2',
    'reference_fasta'
]
def missing = required.findAll { params[it] == null }
if (missing) {
    exit 1, "Missing required param(s): ${missing.join(', ')}. See docs/PHASE1_NOTES.md for a full run example."
}

workflow {

    // Sample IDs match the naming already used throughout docs/ and the run manifest
    // (COLO829 = tumour, COLO829BL = matched normal) -- keep these consistent, Phase 2's
    // Mutect2 tumour/normal pairing and the run_manifest.json samples block both key off them.
    samples_ch = Channel.of(
        tuple('COLO829',   file(params.tumour_reads_1), file(params.tumour_reads_2)),
        tuple('COLO829BL', file(params.normal_reads_1), file(params.normal_reads_2))
    )

    reference_fasta = file(params.reference_fasta)

    SOMATIC(samples_ch, reference_fasta)
}

workflow.onComplete {
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Execution status: ${workflow.success ? 'OK' : 'FAILED'}"
    log.info "Duplicate-marked BAMs: ${params.outdir}/alignment/<sample_id>/<sample_id>.dedup.bam"
    log.info "MultiQC report: ${params.outdir}/qc/multiqc/multiqc_report.html"
}
