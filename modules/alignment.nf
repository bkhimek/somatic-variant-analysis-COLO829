// Module 2 — Alignment
// BWA-MEM2 -> samtools sort -> GATK MarkDuplicates, per plan §6.
// Runs once per sample (tumour, normal each call this independently — see workflows/somatic.nf).
//
// Split into four processes rather than one long pipe: BWA-MEM2 and samtools ship as
// separate single-tool biocontainers (no combined image is used here, to stay consistent
// with the rest of this pipeline's one-container-per-tool convention — see docs/data_sources.md
// for why a "mulled" combined container was deliberately not used). bwa-mem2 mem writes SAM
// to a file; SAMTOOLS_SORT reads that file in a separate process/container. This costs one
// extra intermediate file on disk per sample but keeps every container single-purpose and
// independently upgradable.

process BWA_MEM2_INDEX {
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0'
    publishDir "${params.outdir}/reference/bwa_mem2_index", mode: 'copy'

    input:
    path(reference_fasta)

    output:
    path("${reference_fasta}.*"), emit: index

    script:
    """
    bwa-mem2 index ${reference_fasta}
    """
}

process BWA_MEM2_ALIGN {
    tag "${sample_id}"
    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0'

    input:
    tuple val(sample_id), path(reads1), path(reads2)
    path(reference_fasta)
    path(index_files)

    output:
    tuple val(sample_id), path("${sample_id}.sam"), emit: sam

    script:
    // Read group is required by GATK downstream (Mutect2 in Phase 2 will reject a BAM
    // with no/incomplete @RG). PL:ILLUMINA is correct for PRJEB27698's HiSeq X Ten runs
    // (docs/data_sources.md §1) -- revisit if a different platform's FASTQs are ever used.
    """
    bwa-mem2 mem -t ${task.cpus} \\
        -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}' \\
        ${reference_fasta} ${reads1} ${reads2} > ${sample_id}.sam
    """
}

process SAMTOOLS_SORT {
    tag "${sample_id}"
    // TODO: exact build-hash tag suffix not confirmed from this sandbox (quay.io's tag
    // list is blocked by robots.txt for automated fetching, same limitation noted for
    // CNVkit/hap.py in docs/data_sources.md). Confirm the real tag with
    // `docker pull quay.io/biocontainers/samtools:1.21` on your machine (or browse
    // https://quay.io/repository/samtools/samtools?tab=tags directly) before first run,
    // and update this line to the exact resolved tag.
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    publishDir "${params.outdir}/alignment/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(sam_file)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam ${sam_file}
    samtools index ${sample_id}.sorted.bam
    """
}

process MARK_DUPLICATES {
    tag "${sample_id}"
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.outdir}/alignment/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam), path(sorted_bai)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bam.bai"), emit: bam
    path("${sample_id}.dedup.metrics.txt"), emit: metrics

    script:
    // BQSR is deliberately not run here — plan's documented decision, see
    // docs/data_sources.md §8. Downstream (Phase 2+) reads this dedup BAM directly.
    """
    gatk MarkDuplicates \\
        --INPUT ${sorted_bam} \\
        --OUTPUT ${sample_id}.dedup.bam \\
        --METRICS_FILE ${sample_id}.dedup.metrics.txt \\
        --CREATE_INDEX false

    # GATK MarkDuplicates' own indexer names the file <prefix>.bai, not <prefix>.bam.bai;
    # index with samtools-compatible naming instead so downstream processes (and IGV) find
    # it predictably.
    gatk BuildBamIndex --INPUT ${sample_id}.dedup.bam --OUTPUT ${sample_id}.dedup.bam.bai
    """
}
