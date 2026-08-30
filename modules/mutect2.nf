// Module 4 -- Somatic variant calling (GATK Best Practices: Mutect2 -> LearnReadOrientationModel
// -> FilterMutectCalls). Runs once overall (tumour + normal together in one Mutect2 invocation,
// per GATK's tumour-with-matched-normal design), not per-sample like Modules 1-2.
//
// BQSR is deliberately not applied before this (docs/data_sources.md §8, carried over from the
// plan) -- Mutect2 reads the dedup BAMs from Module 2 directly.

process MUTECT2 {
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.outdir}/mutect2", mode: 'copy'

    input:
    tuple val(tumour_id), path(tumour_bam), path(tumour_bai)
    tuple val(normal_id), path(normal_bam), path(normal_bai)
    path(reference_fasta)
    path(reference_fai)
    path(reference_dict)
    tuple path(panel_of_normals), path(panel_of_normals_index)
    tuple path(germline_resource), path(germline_resource_index)
    path(interval_file)

    output:
    path("unfiltered.vcf.gz"), emit: vcf
    path("unfiltered.vcf.gz.tbi"), emit: vcf_index
    path("unfiltered.vcf.gz.stats"), emit: stats
    path("f1r2.tar.gz"), emit: f1r2

    script:
    // NO_FILE sentinel pattern (see assets/NO_FILE / modules/contamination.nf) -- omitted
    // entirely on the `full` profile, where params.interval_list is null.
    def interval_arg = interval_file.name != 'NO_FILE' ? "-L ${interval_file}" : ''
    // -normal takes the SAMPLE NAME as recorded in the BAM's @RG SM tag, not a file path --
    // modules/alignment.nf's BWA_MEM2_ALIGN sets SM:${sample_id}, so normal_id (COLO829BL)
    // matches exactly what Mutect2 expects to find inside normal_bam's read group.
    //
    // --java-options "-Xmx6g" -- found via execution 2026-08-30: gatk's own default JVM
    // heap-sizing under Docker does not reliably scale to the container's actual memory
    // allocation. The crash log showed `Runtime.totalMemory()=2147483648` (exactly 2GiB
    // self-allocated) inside this process's 8GB container, causing a genuine
    // `java.lang.OutOfMemoryError: Java heap space` while Mutect2 was still just loading the
    // reference sequence dictionary -- not doing any heavy calling work yet. 6g leaves ~2GB of
    // headroom below the 8GB container allocation for off-heap/native memory, JVM overhead, and
    // HTSJDK buffers. This is a systemic property of the gatk wrapper under Docker, not
    // something specific to Mutect2's workload, so the same flag is applied to every GATK
    // invocation in this pipeline (see LEARN_READ_ORIENTATION_MODEL and FILTER_MUTECT_CALLS
    // below, plus modules/contamination.nf and modules/reference_prep.nf).
    """
    gatk --java-options "-Xmx6g" Mutect2 \\
        -R ${reference_fasta} \\
        -I ${tumour_bam} \\
        -I ${normal_bam} \\
        -normal ${normal_id} \\
        --panel-of-normals ${panel_of_normals} \\
        --germline-resource ${germline_resource} \\
        --f1r2-tar-gz f1r2.tar.gz \\
        ${interval_arg} \\
        -O unfiltered.vcf.gz
    """
}

process LEARN_READ_ORIENTATION_MODEL {
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.outdir}/mutect2", mode: 'copy'

    input:
    path(f1r2_tar_gz)

    output:
    path("read-orientation-model.tar.gz"), emit: model

    script:
    // --java-options "-Xmx6g" -- see MUTECT2's script block above for why this is applied to
    // every GATK invocation in this pipeline, not just Mutect2 itself.
    """
    gatk --java-options "-Xmx6g" LearnReadOrientationModel -I ${f1r2_tar_gz} -O read-orientation-model.tar.gz
    """
}

process FILTER_MUTECT_CALLS {
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
    // for MUTECT2 to fully finish writing it, not because it's passed as its own CLI flag).
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
