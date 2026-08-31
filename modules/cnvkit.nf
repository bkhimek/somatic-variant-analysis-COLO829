// Module 6 -- Copy-number calling (Phase 4, added 2026-08-31). Runs CNVkit's whole-genome
// (`-m wgs`) tumour/normal workflow against MARK_DUPLICATES' dedup BAMs, then converts the
// resulting segments into integer copy-number calls.
//
// Module numbering note: earlier docs (docs/data_sources.md, docs/PHASE0_FINDINGS.md) referred
// to this as "Module 7" in a couple of places (with "Module 6" reserved for a VEP/Funcotator
// annotation step) while other places just said "Phase 4 (CNVkit, ...)". The project plan itself
// isn't in this repo, and going back to check it wasn't done for this build -- the user
// explicitly asked for this to be designed from CNVkit's standard workflow rather than blocked
// on retrieving the plan's exact wording (2026-08-31). Called "Module 6" here as a working label
// only; renumber freely once/if the plan's real numbering is confirmed -- it has no effect on
// how the code runs.
//
// No separate "PREPARE_..." auto-build step here, unlike modules/reference_prep.nf's .fai/.dict
// pattern -- CNVkit's own `batch` command builds its per-run reference (from the normal sample)
// and access BED (from the reference FASTA's N-gaps) internally, on the fly, each run. Unlike
// the bwa-mem2 index or GATK .fai/.dict, these aren't expensive, shared, reusable-across-runs
// artifacts worth caching outside the process -- so there's nothing to auto-detect-or-skip here.
//
// Known, deliberate scope limits for this first cut (see docs/PHASE4_NOTES.md for the full
// reasoning -- same "prove the wiring first" split established in Phases 1-3):
//   - No `--annotate refFlat.txt`. CNVkit's own docs call a gene-annotation database necessary
//     for *readable per-target gene labels* in the output, not for the command to run at all --
//     bins will carry generic/positional labels instead of real gene symbols. Sourcing and
//     verifying a GRCh38 refFlat.txt is deferred rather than rushed in for a first cut whose
//     query BAM has almost no real coverage to label anyway.
//   - No `--scatter`/`--diagram`. Both are valid `batch` flags and would be nice-to-have QC
//     plots, but this project just got burned (Phase 3, som.py's `--happy-stats`) by adding an
//     extra output-producing flag without checking what it actually needs to succeed -- not
//     repeating that here for a first cut. Can be added once the base call path is proven to run
//     cleanly against real data.
//   - `-y`/`--male-reference` IS included, not deferred: verified (not guessed) that the COLO829
//     donor is male, 45 years old, per ATCC's own COLO 829 (CRL-1974) product page -- getting
//     chrX/chrY ploidy handling wrong from the start would silently bias every sex-chromosome
//     call, so this one was worth confirming before writing any code rather than defaulting to
//     CNVkit's female-reference assumption.
//   - `--drop-low-coverage` IS included: CNVkit's own docs recommend it specifically "for tumor
//     samples" (drops bins at/near zero read depth rather than segmenting noise), which matches
//     our situation more than usual -- the only tumour BAM that exists so far (Phase 2's
//     `dev`-profile result) comes from a 10,000-read-pair subsample spread across the whole
//     genome, so most whole-genome bins will have effectively zero coverage.
//
// First-run finding (2026-08-31, see docs/PHASE4_NOTES.md "First run -- findings" for the full
// trail): the first real execution failed with "Missing output file(s) `*.cns`" -- CNVKIT_BATCH
// exited 0 but never wrote a .cns, because `fix` logged "Keeping 0 of 58496 bins" and wrote an
// empty (0-region) .cnr for `segment` to work from. Traced (via CNVkit's own `fix.py` source, not
// guessed) to `mask_bad_bins()` -- an always-on filter against the REFERENCE's own per-bin
// log2/spread QC, unrelated to `--drop-low-coverage` (that flag only affects `segment`'s handling
// of the tumour side, and never got the chance to matter here). The reference itself, built from
// the equally-sparse normal BAM, had already logged "100.0% bins failed filters" during
// `reference` construction -- with CNVkit's autobin-computed ~53kb bins and this subsample's
// ~20,000 reads spread across the whole 3.1Gb genome (~0.001x effective coverage), the expected
// reads per bin is a fraction of one, so essentially every bin's log2/spread statistics are
// degenerate. This is a genuine statistical floor, not a flag or container bug -- the same shape
// as Phase 1's "duplication rate needs real depth" and Phase 2's "unsharded genome-wide Mutect2
// doesn't fit this machine" findings. Fix: `params.cnvkit_target_avg_size` (set to 10Mb on the
// `dev` profile only, see nextflow.config) forces far larger, far fewer bins, concentrating the
// same read count enough to plausibly clear the reference-quality filter. Left null on `full` --
// real WGS depth shouldn't need this override, since CNVkit's own autobin is designed for it.
//
// Second-run finding (2026-08-31, see docs/PHASE4_NOTES.md "Second run -- findings" for the full
// trail): with the bin-size fix above, CNVKIT_BATCH itself succeeded -- but CNVKIT_CALL then
// failed, because CNVKIT_BATCH's work directory actually contained THREE `*.cns`-suffixed files
// (`<prefix>.cns`, `<prefix>.call.cns`, `<prefix>.bintest.cns`), not the one plain segmented file
// this module assumed. The broad `path("*.cns")` glob matched all three, handing CNVKIT_CALL a
// 3-element list instead of one file, which Nextflow interpolated into the command line as three
// space-separated positional arguments -- `cnvkit.py call` rejected the extras as "unrecognized
// arguments." Checked against CNVkit's own current (master-branch) `commands.py` and
// `segmentation/__init__.py` source: neither shows `batch` invoking `call`/`bintest`
// automatically, so this may be specific to the older pinned 0.9.10 release (bioconda's cnvkit
// recipe is already at 0.9.13 as of this writing) or some other version/flag interaction not
// visible from current source -- left as a confirmed-by-execution fact rather than a fully
// root-caused one. Fix: made the `cnr`/`cns` outputs exact filenames (CNVkit's own naming
// convention from the input BAM's basename, already known and used correctly elsewhere in this
// pipeline) instead of glob patterns, so CNVKIT_CALL only ever receives the one plain segmented
// file regardless of what else `batch` happens to produce alongside it. Also split
// CNVKIT_BATCH/CNVKIT_CALL into separate `publishDir` subfolders (`cnvkit/batch`, `cnvkit/call`)
// -- otherwise CNVKIT_BATCH's own mystery `<prefix>.call.cns` and CNVKIT_CALL's genuine
// `<prefix>.call.cns` would silently collide in the same output directory.

process CNVKIT_BATCH {
    tag { tumour_id }
    // Container tag inherited from Phase 0 research (docs/data_sources.md §7): the plan's
    // original choice. Unlike htslib and hap.py in Phase 3 (both guessed wrong, since quay.io
    // blocks automated tag-list fetching from this sandbox), this one was confirmed correct by
    // a real execution 2026-08-31 -- pulled and ran without any manifest-not-found error.
    container 'quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0'
    publishDir "${params.outdir}/cnvkit/batch", mode: 'copy'

    input:
    tuple val(tumour_id), path(tumour_bam), path(tumour_bai)
    tuple val(normal_id), path(normal_bam), path(normal_bai)
    path(reference_fasta)
    path(reference_fai)
    val(target_avg_size)   // params.cnvkit_target_avg_size -- null on `full` (let autobin decide), a large fixed value on `dev` (see nextflow.config and this module's header comment for why)

    output:
    path("${tumour_bam.baseName}.cnr"), emit: cnr
    path("${tumour_bam.baseName}.cns"), emit: cns   // exact filename, not a glob -- see "Second-run finding" above for why
    path("*.cnn"), emit: coverage_files, optional: true
    path("*"), emit: all_outputs   // catch-all for reference.cnn plus whatever else batch produces (e.g. the mystery .call.cns/.bintest.cns from the finding above) -- kept for reference, not consumed downstream

    script:
    // -m wgs: whole-genome mode -- treats the reference's sequencing-accessible regions as
    // CNVkit's "targets" (auto-computed on the fly, confirmed via CNVkit's own docs, since no
    // -t/--targets or -g/--access file is supplied) rather than a captured-exome/panel BED.
    // -y: see module header comment -- COLO829's donor is male, confirmed via ATCC CRL-1974.
    // --drop-low-coverage: see module header comment -- recommended by CNVkit for tumor samples,
    // and especially relevant given how sparse the only tumour BAM so far actually is.
    // --target-avg-size: only added when target_avg_size is set (dev profile) -- see the
    // "First-run finding" in this module's header comment for why the dev profile needs this and
    // the full profile deliberately doesn't.
    // -p ${task.cpus}: parallelise CNVkit's own per-chromosome/per-sample steps, same convention
    // as bwa-mem2's -t and Nextflow's own task.cpus elsewhere in this pipeline.
    // -d .: write outputs into the process's own work directory (Nextflow's usual convention --
    // publishDir above copies them out), rather than CNVkit's own default of the current
    // directory (which would be the same thing here, but this makes the choice explicit).
    def target_avg_size_arg = target_avg_size ? "--target-avg-size ${target_avg_size}" : ''
    """
    cnvkit.py batch ${tumour_bam} \\
        --normal ${normal_bam} \\
        -m wgs \\
        -f ${reference_fasta} \\
        -y \\
        --drop-low-coverage \\
        ${target_avg_size_arg} \\
        -p ${task.cpus} \\
        -d .
    """
}

process CNVKIT_CALL {
    tag { cns_file.baseName }
    container 'quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0'
    publishDir "${params.outdir}/cnvkit/call", mode: 'copy'

    input:
    path(cns_file)   // CNVKIT_BATCH.out.cns -- now an exact single file, not a glob (see that process's output block)

    output:
    path("*.call.cns"), emit: call_cns

    script:
    // -y here too: CNVkit's own docs are explicit that if -y/--male-reference was used to build
    // the reference (as it was in CNVKIT_BATCH above), the same flag must be repeated on `call`
    // for consistent sex-chromosome ploidy handling -- easy to silently get wrong by omission,
    // so called out here rather than assumed carried-over.
    """
    cnvkit.py call ${cns_file} -y -o ${cns_file.baseName}.call.cns
    """
}
