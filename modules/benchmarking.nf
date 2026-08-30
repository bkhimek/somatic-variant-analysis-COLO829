// Module 5 -- Benchmarking (Phase 3, added 2026-08-30). Compares FILTER_MUTECT_CALLS' output
// against the NYGC COLO-829/COLO-829BL open somatic VCF (decided as the primary SNV/indel truth
// set in Phase 0 -- docs/PHASE0_FINDINGS.md, docs/data_sources.md) using hap.py, the standard
// GA4GH-recommended tool for benchmarking a variant caller's precision/recall against a truth
// set.
//
// Known, deliberate scope limits for this first cut (see docs/PHASE3_NOTES.md for the full
// reasoning):
//   - No confidence-region BED restriction -- NYGC doesn't publish one for SNV/indel (only CNV
//     bed + SV bedpe exist per cell line, confirmed in Phase 0), so hap.py runs unrestricted.
//   - Default hap.py comparison engine (xcmp), not vcfeval -- avoids pulling in RTG Tools/an SDF
//     reference just for this first pass; xcmp is hap.py's own default and needs nothing beyond
//     what's already in the container.
//   - No --pass-only override -- hap.py's own default already only counts PASS variants in both
//     truth and query, which matches what we want (FILTER_MUTECT_CALLS' filtered.vcf.gz marks
//     non-PASS calls with real filter reasons; those should be excluded from being counted as
//     positive calls, not force-included).

process PREPARE_TRUTH_VCF {
    // Found via execution 2026-08-30: quay.io/biocontainers/htslib:1.21--h566b1c6_0 doesn't
    // exist ("manifest unknown") -- that tag was guessed by pattern-matching the samtools image
    // tag's build-hash convention rather than actually checked. Verified via research afterwards
    // rather than guessing again: reusing quay.io/biocontainers/samtools:1.21--h50ea8bc_0 (already
    // confirmed working in this pipeline since Phase 1) is both simpler and correct -- its
    // bioconda recipe pulls in htslib as a real runtime dependency (via htslib's own
    // `run_exports`, not just a build-time link), which installs the bgzip/tabix binaries into
    // the same conda environment alongside samtools itself. Avoids pulling a second, separately
    // unverified image for what both tools already ship together. See docs/PHASE3_NOTES.md.
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    publishDir "${params.outdir}/benchmarking/truth_set", mode: 'copy'

    input:
    path(truth_vcf_plain)

    output:
    tuple path("${truth_vcf_plain.baseName}.vcf.gz"), path("${truth_vcf_plain.baseName}.vcf.gz.tbi"), emit: truth_vcf_indexed

    script:
    // NYGC ships the truth set as a plain (non-bgzip) uncompressed VCF (confirmed in Phase 0,
    // docs/data_sources.md sec 3) -- hap.py wants bgzip+tabix like every other VCF input in this
    // pipeline, so this mirrors modules/reference_prep.nf's INDEX_VCF pattern rather than
    // reinventing it, just with an added bgzip step since this file isn't gzipped at all yet
    // (the resource VCFs INDEX_VCF handles are already bgzipped, just missing their .tbi).
    """
    bgzip -c ${truth_vcf_plain} > ${truth_vcf_plain.baseName}.vcf.gz
    tabix -p vcf ${truth_vcf_plain.baseName}.vcf.gz
    """
}

process HAPPY_BENCHMARK {
    container 'quay.io/biocontainers/hap.py:0.3.15-0'
    publishDir "${params.outdir}/benchmarking", mode: 'copy'

    input:
    tuple path(truth_vcf), path(truth_vcf_index)
    path(query_vcf)
    path(query_vcf_index)
    path(reference_fasta)
    path(reference_fai)

    output:
    path("happy.summary.csv"), emit: summary
    path("happy.extended.csv"), emit: extended
    path("happy.*"), emit: all_outputs

    script:
    // hap.py's own defaults do the right thing here (see module header comment): PASS-only
    // counting on both sides, xcmp comparison engine, no confidence-region restriction since
    // none exists for this truth set. `-o happy` sets the output filename prefix -- hap.py
    // writes happy.summary.csv (the headline precision/recall/F1 table), happy.extended.csv
    // (broken down further), happy.vcf.gz (annotated combined callset), and a few others.
    """
    hap.py \\
        ${truth_vcf} \\
        ${query_vcf} \\
        -r ${reference_fasta} \\
        -o happy
    """
}
