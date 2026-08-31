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

process CNVKIT_BATCH {
    tag { tumour_id }
    // Container tag inherited from Phase 0 research (docs/data_sources.md §7): the plan's
    // original choice, and the bioconda recipe for cnvkit 0.9.10 build 0 was confirmed to exist
    // -- but the exact quay.io tag string itself was NOT independently verified, because quay.io
    // blocks automated tag-list fetching from this sandbox (the same limitation that already
    // produced two wrong-tag bugs this project has hit and fixed by execution -- htslib in
    // Phase 3's "First run -- findings", hap.py in its "Second run -- findings"). Given that
    // history, don't run this without confirming the tag yourself first:
    //     docker pull quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0
    // See docs/PHASE4_NOTES.md for what to do if that pull fails.
    container 'quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0'
    publishDir "${params.outdir}/cnvkit", mode: 'copy'

    input:
    tuple val(tumour_id), path(tumour_bam), path(tumour_bai)
    tuple val(normal_id), path(normal_bam), path(normal_bai)
    path(reference_fasta)
    path(reference_fai)

    output:
    path("*.cnr"), emit: cnr
    path("*.cns"), emit: cns
    path("*.cnn"), emit: coverage_files, optional: true
    path("*"), emit: all_outputs

    script:
    // -m wgs: whole-genome mode -- treats the reference's sequencing-accessible regions as
    // CNVkit's "targets" (auto-computed on the fly, confirmed via CNVkit's own docs, since no
    // -t/--targets or -g/--access file is supplied) rather than a captured-exome/panel BED.
    // -y: see module header comment -- COLO829's donor is male, confirmed via ATCC CRL-1974.
    // --drop-low-coverage: see module header comment -- recommended by CNVkit for tumor samples,
    // and especially relevant given how sparse the only tumour BAM so far actually is.
    // -p ${task.cpus}: parallelise CNVkit's own per-chromosome/per-sample steps, same convention
    // as bwa-mem2's -t and Nextflow's own task.cpus elsewhere in this pipeline.
    // -d .: write outputs into the process's own work directory (Nextflow's usual convention --
    // publishDir above copies them out), rather than CNVkit's own default of the current
    // directory (which would be the same thing here, but this makes the choice explicit).
    """
    cnvkit.py batch ${tumour_bam} \\
        --normal ${normal_bam} \\
        -m wgs \\
        -f ${reference_fasta} \\
        -y \\
        --drop-low-coverage \\
        -p ${task.cpus} \\
        -d .
    """
}

process CNVKIT_CALL {
    tag { cns_file.baseName }
    container 'quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0'
    publishDir "${params.outdir}/cnvkit", mode: 'copy'

    input:
    path(cns_file)

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
