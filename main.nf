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

workflow {

    // ---- Required params: fail fast with a clear message rather than a cryptic NPE mid-DAG ----
    // Moved inside the workflow block 2026-08-30 -- Nextflow 26.x's stricter parser rejects
    // top-level imperative statements (def/if) mixed with declarative script elements
    // (include/workflow), the same category of break already hit and fixed in nextflow.config.
    // See docs/PHASE1_NOTES.md for detail.
    def required = [
        'tumour_reads_1', 'tumour_reads_2',
        'normal_reads_1', 'normal_reads_2',
        'reference_fasta'
    ]
    def missing = required.findAll { params[it] == null }
    if (missing) {
        exit 1, "Missing required param(s): ${missing.join(', ')}. See docs/PHASE1_NOTES.md for a full run example."
    }

    // Sample IDs match the naming already used throughout docs/ and the run manifest
    // (COLO829 = tumour, COLO829BL = matched normal) -- keep these consistent, Phase 2's
    // Mutect2 tumour/normal pairing and the run_manifest.json samples block both key off them.
    samples_ch = Channel.of(
        tuple('COLO829',   file(params.tumour_reads_1), file(params.tumour_reads_2)),
        tuple('COLO829BL', file(params.normal_reads_1), file(params.normal_reads_2))
    )

    reference_fasta = file(params.reference_fasta)

    SOMATIC(samples_ch, reference_fasta)

    // A custom workflow.onComplete{} summary block used to live here (added, moved inside this
    // block for Nextflow 26.x's parser, then wrapped in try/catch -- see git history / earlier
    // revisions of docs/PHASE1_NOTES.md for the full saga). Removed entirely 2026-08-30: it kept
    // throwing "ERROR ~ Failed to invoke `workflow.onComplete` event handler" even wrapped in
    // try/catch, which means the failure happens in Nextflow's own mechanism for invoking a
    // handler registered from inside the entry workflow's execution scope -- not in anything our
    // closure's body did (a try/catch around simple log.info calls can't be the culprit; the
    // wrapper wasn't even reaching our log.warn fallback). Since the run's actual outcome is
    // already reported by Nextflow's own built-in completion summary (Completed at / Duration /
    // Succeeded / Cached -- printed automatically, no config needed) and this block only ever
    // added a few redundant path reminders, it wasn't worth chasing further. Output paths are
    // documented in "How to run it" below instead, since they're static regardless of run.
}
