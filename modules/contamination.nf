// Module 3 -- Contamination estimation (GATK Best Practices: GetPileupSummaries + CalculateContamination).
//
// Unlike Modules 1-2, which process tumour and normal identically and independently,
// CalculateContamination genuinely needs BOTH samples' pileup tables together (it estimates
// contamination in the tumour using the matched normal as a reference) -- see workflows/somatic.nf
// for how the two samples are branched apart, fed through the same GET_PILEUP_SUMMARIES process,
// then rejoined for CALCULATE_CONTAMINATION.

process GET_PILEUP_SUMMARIES {
    tag { sample_id }
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir { "${params.outdir}/contamination/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path(reference_fasta)
    path(reference_fai)
    path(reference_dict)
    tuple path(common_biallelic_sites), path(common_biallelic_sites_index)
    path(interval_file)

    output:
    tuple val(sample_id), path("${sample_id}.pileups.table"), emit: pileups

    script:
    // Found via execution 2026-08-30: unlike Mutect2, GATK's GetPileupSummaries genuinely
    // REQUIRES -L/--intervals -- omitting it entirely (the original NO_FILE-sentinel design,
    // matching how the same sentinel makes -L optional for Mutect2 in modules/mutect2.nf) fails
    // with "A USER ERROR has occurred: Argument intervals was missing: Argument 'intervals' is
    // required". Standard GATK Best Practices fix: always restrict to the common-sites VCF
    // itself as the baseline interval set (those are the only positions GetPileupSummaries ever
    // needs to check), and additionally intersect with a real interval_list when one is given
    // (the dev-profile driver-gene BED) via --interval-set-rule INTERSECTION, narrowing further
    // to just the genes' overlap with common SNP sites rather than replacing the sites
    // restriction outright.
    def interval_args = interval_file.name != 'NO_FILE'
        ? "-L ${interval_file} -L ${common_biallelic_sites} --interval-set-rule INTERSECTION"
        : "-L ${common_biallelic_sites}"
    """
    gatk GetPileupSummaries \\
        -I ${bam} \\
        -R ${reference_fasta} \\
        -V ${common_biallelic_sites} \\
        ${interval_args} \\
        -O ${sample_id}.pileups.table
    """
}

process CALCULATE_CONTAMINATION {
    container 'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.outdir}/contamination", mode: 'copy'

    input:
    path(tumour_pileups)
    path(normal_pileups)

    output:
    path("contamination.table"), emit: contamination_table
    path("segments.table"), emit: segments_table

    script:
    // Threshold check logged (not enforced -- exit code stays 0 regardless) per
    // conf/resources.config's explicit note that these defaults aren't final until inspected
    // against real COLO829 output. This just surfaces the number prominently in the log
    // rather than requiring you to go open contamination.table by hand.
    """
    gatk CalculateContamination \\
        -I ${tumour_pileups} \\
        -matched ${normal_pileups} \\
        -O contamination.table \\
        --tumor-segmentation segments.table

    CONTAM_FRACTION=\$(awk 'NR==2 {print \$2}' contamination.table)
    WARN_THRESHOLD=${params.qc.contamination.warn_threshold}
    FAIL_THRESHOLD=${params.qc.contamination.fail_threshold}
    echo "Contamination fraction: \${CONTAM_FRACTION} (warn > \${WARN_THRESHOLD}, fail > \${FAIL_THRESHOLD})"
    awk -v c="\${CONTAM_FRACTION}" -v warn="\${WARN_THRESHOLD}" -v fail="\${FAIL_THRESHOLD}" 'BEGIN {
        if (c > fail)      { print "FAIL: contamination fraction " c " exceeds fail threshold " fail }
        else if (c > warn) { print "WARN: contamination fraction " c " exceeds warn threshold " warn }
        else                { print "OK: contamination fraction " c " is within threshold" }
    }'
    """
}
