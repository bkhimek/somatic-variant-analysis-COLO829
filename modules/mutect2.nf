// Module 4 -- Somatic variant calling (GATK Best Practices: Mutect2 -> LearnReadOrientationModel
// -> FilterMutectCalls). Runs once overall (tumour + normal together in one Mutect2 invocation,
// per GATK's tumour-with-matched-normal design), not per-sample like Modules 1-2.
//
// BQSR is deliberately not applied before this (docs/data_sources.md §8, carried over from the
// plan) -- Mutect2 reads the dedup BAMs from Module 2 directly.
//
// Interval scatter/gather (added 2026-09-01, see docs/MUTECT2_SCATTERGATHER_NOTES.md for the
// full research trail): Phase 2 (docs/PHASE2_NOTES.md) found that a single unsharded Mutect2
// invocation against the whole genome doesn't fit this machine at any memory size tried, and
// deferred a real fix until full-genome data actually made it necessary. That point has now
// arrived. GATK's own docs, Broad's production `mutect2.wdl`, and nf-core/sarek's actual source
// all confirm the same shape of fix: split the calling regions into N interval-list shards
// (`SplitIntervals`), run Mutect2 once per shard, then merge the per-shard VCFs (`MergeVcfs`,
// NOT `GatherVcfs` -- Broad's own WDL explicitly avoids `GatherVcfs` because it doesn't reliably
// support index creation on block-compressed VCFs) and `.stats` files (`MergeMutectStats`)
// before `FilterMutectCalls` ever sees them. Per-shard `f1r2.tar.gz` files are NOT pre-merged --
// `LearnReadOrientationModel` natively accepts multiple `-I` inputs in one call and does the
// combining itself, exactly as Broad's WDL and sarek both do it.
//
// MUTECT2 itself is now a scattered process (one task per interval shard, wired in
// workflows/somatic.nf) instead of a single call across the whole/dev-restricted genome. Its
// outputs use exact, shard-derived filenames rather than a generic name or a glob -- applying
// the lesson from Phase 4's CNVKIT_CALL bug (docs/PHASE4_NOTES.md's "Second run -- findings":
// a `path("*.cns")` glob matching more files than expected silently broke a downstream single-
// file consumer) proactively here, rather than waiting to hit the same mistake again.

process SCATTER_INTERVALS_BY_NS {
    // Only invoked on the `full` profile (see workflows/somatic.nf) -- when there's no
    // interval_list restriction, the whole genome needs an N-gap-aware ACGT-only interval list
    // before it's safe to split by SplitIntervals. Splitting blindly by raw base count (which
    // SplitIntervals does if handed no -L at all) can cut a shard boundary through an assembly
    // gap; this is exactly how Broad's own public `wgs_calling_regions.hg38.interval_list` is
    // built (docs/MUTECT2_SCATTERGATHER_NOTES.md §1). On `dev`, params.interval_list is already
    // a small curated real-gene BED -- no N-gap scan needed there, so this process is skipped
    // entirely (workflows/somatic.nf branches on params.interval_list, not on profile name).
    container 'broadinstitute/gatk:4.5.0.0'

    input:
    path(reference_fasta)
    path(reference_fai)
    path(reference_dict)

    output:
    path("acgt.interval_list"), emit: acgt_intervals

    script:
    // -OT ACGT: only the sequenceable (non-N) blocks -- confirmed via GATK's own
    // ScatterIntervalsByNs doc (docs/MUTECT2_SCATTERGATHER_NOTES.md §1), default is BOTH
    // (N-blocks and ACGT-blocks together), which we don't want here.
    """
    gatk --java-options "-Xmx6g" ScatterIntervalsByNs \\
        -R ${reference_fasta} \\
        -O acgt.interval_list \\
        -OT ACGT
    """
}

process SPLIT_INTERVALS {
    container 'broadinstitute/gatk:4.5.0.0'

    input:
    path(reference_fasta)
    path(reference_fai)
    path(reference_dict)
    path(calling_intervals)
    val(scatter_count)

    output:
    path("split_intervals/*-scattered.interval_list"), emit: shards

    script:
    // --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION: keeps each input interval
    // intact rather than splitting inside one (GATK's own SplitIntervals doc recommends this
    // for assembly-based callers like Mutect2, to avoid analytical edge artifacts at a shard
    // boundary that falls mid-interval -- docs/MUTECT2_SCATTERGATHER_NOTES.md §1).
    """
    mkdir -p split_intervals
    gatk --java-options "-Xmx6g" SplitIntervals \\
        -R ${reference_fasta} \\
        -L ${calling_intervals} \\
        --scatter-count ${scatter_count} \\
        --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \\
        -O split_intervals
    """
}

process MUTECT2 {
    tag { interval_shard.baseName }
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.outdir}/mutect2/scatter", mode: 'copy'
    // cpus dropped to 1 (from the flat 4-cpu profile default) -- GATK's own Mutect2 FAQ says
    // one CPU per instance is enough, since parallelism here is meant to come from running many
    // shards concurrently, not from multithreading a single one (docs/MUTECT2_SCATTERGATHER_
    // NOTES.md §5). This lets more shards actually run at once within the executor's shared
    // cpu/memory pool (nextflow.config) instead of one 4-cpu task blocking the rest.
    cpus 1

    input:
    tuple val(tumour_id), path(tumour_bam), path(tumour_bai)
    tuple val(normal_id), path(normal_bam), path(normal_bai)
    path(reference_fasta)
    path(reference_fai)
    path(reference_dict)
    tuple path(panel_of_normals), path(panel_of_normals_index)
    tuple path(germline_resource), path(germline_resource_index)
    path(interval_shard)

    output:
    path("${interval_shard.baseName}.unfiltered.vcf.gz"), emit: vcf
    path("${interval_shard.baseName}.unfiltered.vcf.gz.tbi"), emit: vcf_index
    path("${interval_shard.baseName}.unfiltered.vcf.gz.stats"), emit: stats
    path("${interval_shard.baseName}.f1r2.tar.gz"), emit: f1r2

    script:
    // -normal takes the SAMPLE NAME as recorded in the BAM's @RG SM tag, not a file path --
    // modules/alignment.nf's BWA_MEM2_ALIGN sets SM:${sample_id}, so normal_id (COLO829BL)
    // matches exactly what Mutect2 expects to find inside normal_bam's read group.
    //
    // --java-options "-Xmx6g" -- found via execution 2026-08-30 (see git history/PHASE2_NOTES.md):
    // gatk's own default JVM heap-sizing under Docker does not reliably scale to the container's
    // actual memory allocation. Applied to every GATK invocation in this pipeline for the same
    // reason, unchanged by this scatter/gather rework.
    //
    // Exact, shard-derived output filenames (not a fixed name, not a glob) -- every shard task
    // writes into its own isolated Nextflow work directory, so same-named outputs across shards
    // don't collide there, but publishDir copies every task's output into ONE shared results
    // directory -- without a unique name per shard, N tasks would all try to publish
    // "unfiltered.vcf.gz" and silently overwrite each other (the same category of mistake as
    // Phase 4's CNVKIT_CALL glob bug, just from a different direction -- collision on the way
    // out instead of ambiguity on the way in).
    """
    gatk --java-options "-Xmx6g" Mutect2 \\
        -R ${reference_fasta} \\
        -I ${tumour_bam} \\
        -I ${normal_bam} \\
        -normal ${normal_id} \\
        --panel-of-normals ${panel_of_normals} \\
        --germline-resource ${germline_resource} \\
        --f1r2-tar-gz ${interval_shard.baseName}.f1r2.tar.gz \\
        -L ${interval_shard} \\
        -O ${interval_shard.baseName}.unfiltered.vcf.gz
    """
}

process MERGE_VCFS {
    // MergeVcfs, not GatherVcfs -- Broad's own mutect2.wdl explicitly picks MergeVcfs "so we can
    // create indices" (docs/MUTECT2_SCATTERGATHER_NOTES.md §2 quotes their exact comment: as of
    // that WDL's writing, GatherVcfs index creation wasn't reliably supported for block-
    // compressed VCFs). Same choice made here for the same reason, not independently guessed.
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.outdir}/mutect2", mode: 'copy'

    input:
    path(vcfs)
    // vcf_indices isn't referenced directly in the command below (MergeVcfs finds each .tbi by
    // naming convention next to its .vcf.gz, same auto-discovery FILTER_MUTECT_CALLS already
    // relies on for its own .stats input) -- declaring it as an input just makes Nextflow
    // actually stage every .tbi into this task's work directory alongside its .vcf.gz, the same
    // "explicit dependency, not a CLI flag" pattern already used in FILTER_MUTECT_CALLS below.
    path(vcf_indices)

    output:
    path("unfiltered.vcf.gz"), emit: vcf
    path("unfiltered.vcf.gz.tbi"), emit: vcf_index

    script:
    // -I repeated once per shard VCF -- confirmed via GATK's own MergeVcfs doc
    // (docs/MUTECT2_SCATTERGATHER_NOTES.md §2), not the single -I-with-a-.list-file alternative
    // it also supports; repeating -I is simpler to build from a Nextflow-collected list.
    def inputs_arg = vcfs.collect { "-I ${it}" }.join(' ')
    """
    gatk --java-options "-Xmx6g" MergeVcfs \\
        ${inputs_arg} \\
        -O unfiltered.vcf.gz
    """
}

process MERGE_MUTECT_STATS {
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.outdir}/mutect2", mode: 'copy'

    input:
    path(stats_files)

    output:
    // Named to match MERGE_VCFS' "unfiltered.vcf.gz" exactly (<vcf-name>.stats) -- FILTER_
    // MUTECT_CALLS below relies on GATK's own naming-convention auto-discovery (finds
    // "unfiltered.vcf.gz.stats" sitting next to "unfiltered.vcf.gz" in its own task work
    // directory, since Nextflow stages both declared path inputs into the same directory
    // regardless of which upstream process produced each one) rather than an explicit --stats
    // flag -- same pattern the single-shot version of this pipeline already used.
    path("unfiltered.vcf.gz.stats"), emit: stats

    script:
    // -stats repeated once per shard's .stats file -- confirmed via GATK's own worked example
    // for scattered Mutect2 (docs/MUTECT2_SCATTERGATHER_NOTES.md §2/§4).
    def stats_arg = stats_files.collect { "-stats ${it}" }.join(' ')
    """
    gatk --java-options "-Xmx6g" MergeMutectStats \\
        ${stats_arg} \\
        -O unfiltered.vcf.gz.stats
    """
}

process LEARN_READ_ORIENTATION_MODEL {
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.outdir}/mutect2", mode: 'copy'

    input:
    // Now a collected list of every shard's f1r2.tar.gz, not one file -- LearnReadOrientationModel
    // is NOT run per-shard-then-merged; it natively accepts multiple -I inputs in a single call
    // and does the combining itself. Confirmed via GATK's own community guidance for scattered
    // Mutect2 runs (docs/MUTECT2_SCATTERGATHER_NOTES.md §2/§4): "if you are scattering Mutect2
    // ... you must input the --f1r2-tar-gz output from each Mutect2 scatter to
    // LearnReadOrientationModel" -- this is the same tool call as before this rework, just now
    // given N inputs instead of 1 (N=1 still works identically on `dev` if its scatter count is
    // left at 1, and behaves the same either way since -I repeated once is just -I once).
    path(f1r2_tar_gz_files)

    output:
    path("read-orientation-model.tar.gz"), emit: model

    script:
    // --java-options "-Xmx6g" -- see MUTECT2's script block above for why this is applied to
    // every GATK invocation in this pipeline, not just Mutect2 itself.
    def inputs_arg = f1r2_tar_gz_files.collect { "-I ${it}" }.join(' ')
    """
    gatk --java-options "-Xmx6g" LearnReadOrientationModel ${inputs_arg} -O read-orientation-model.tar.gz
    """
}

process FILTER_MUTECT_CALLS {
    // Unchanged from the pre-scatter version -- it always took one merged-per-sample VCF/stats/
    // model/contamination-table set as input; the only difference now is that MERGE_VCFS/
    // MERGE_MUTECT_STATS actually produce those from N shards instead of MUTECT2 producing them
    // directly from a single unsharded call. Broad's own WDL and nf-core/sarek both gather
    // everything to exactly this shape (one-per-sample) before their own FilterMutectCalls call
    // too (docs/MUTECT2_SCATTERGATHER_NOTES.md §2/§3).
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.outdir}/mutect2", mode: 'copy'

    input:
    path(unfiltered_vcf)
    path(unfiltered_vcf_index)
    path(stats)
    path(reference_fasta)
    path(reference_fai)
    path(reference_dict)
    path(contamination_table)
    path(segments_table)
    path(read_orientation_model)

    output:
    path("filtered.vcf.gz"), emit: vcf
    path("filtered.vcf.gz.tbi"), emit: vcf_index
    path("filtered.vcf.gz.filteringStats.tsv"), emit: filtering_stats

    script:
    // GATK finds unfiltered_vcf.stats automatically alongside -V by naming convention (the
    // `stats` input above just makes the dependency explicit to Nextflow so this process waits
    // for MERGE_MUTECT_STATS to fully finish writing it, not because it's passed as its own CLI
    // flag).
    //
    // --java-options "-Xmx6g" -- see MUTECT2's script block above for why this is applied to
    // every GATK invocation in this pipeline, not just Mutect2 itself.
    """
    gatk --java-options "-Xmx6g" FilterMutectCalls \\
        -R ${reference_fasta} \\
        -V ${unfiltered_vcf} \\
        --contamination-table ${contamination_table} \\
        --tumor-segmentation ${segments_table} \\
        --ob-priors ${read_orientation_model} \\
        -O filtered.vcf.gz
    """
}
