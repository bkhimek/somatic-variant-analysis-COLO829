// Module 5 -- Benchmarking (Phase 3, added 2026-08-30). Compares FILTER_MUTECT_CALLS' output
// against the NYGC COLO-829/COLO-829BL open somatic VCF (decided as the primary SNV/indel truth
// set in Phase 0 -- docs/PHASE0_FINDINGS.md, docs/data_sources.md).
//
// Uses som.py, NOT hap.py, despite both shipping in the same package/container -- found via
// execution 2026-08-30 (see docs/PHASE3_NOTES.md): hap.py is GA4GH's tool for GERMLINE
// small-variant benchmarking and keys off the VCF's GT (genotype) field to decide whether a
// record is "non-reference." NYGC's somatic truth VCF (and Mutect2's own output) carry AD/DP/AF
// per-sample annotation but no GT field at all -- with hap.py, every record silently reads as
// homozygous-reference, producing "Non-reference VCF records: 0" despite 43,192 real truth
// records, and hap.py then hard-fails with "Input files/regions do not contain variants." The
// same Illumina/hap.py package ships `som.py` specifically for this: somatic comparison based on
// ALT-allele presence (via bcftools isec under the hood), no GT/haplotype matching required --
// exactly what a Mutect2-vs-somatic-truth-set comparison needs. No new container was needed;
// som.py is already present in the same hap.py image alongside hap.py itself.
//
// Known, deliberate scope limits for this first cut (see docs/PHASE3_NOTES.md for the full
// reasoning):
//   - No confidence-region BED restriction -- NYGC doesn't publish one for SNV/indel (only CNV
//     bed + SV bedpe exist per cell line, confirmed in Phase 0), so som.py runs unrestricted.
//   - `-N` (normalize both VCFs before comparing) -- som.py's own documented recommendation for
//     a first run.
//   - `--happy-stats` -- adds a friendlier som.summary.csv table alongside som.py's raw
//     som.stats.csv, giving the same kind of headline precision/recall view hap.py's
//     summary.csv would have provided.

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
    // docs/data_sources.md sec 3) -- som.py wants bgzip+tabix like every other VCF input in this
    // pipeline, so this mirrors modules/reference_prep.nf's INDEX_VCF pattern rather than
    // reinventing it, just with an added bgzip step since this file isn't gzipped at all yet
    // (the resource VCFs INDEX_VCF handles are already bgzipped, just missing their .tbi).
    """
    bgzip -c ${truth_vcf_plain} > ${truth_vcf_plain.baseName}.vcf.gz
    tabix -p vcf ${truth_vcf_plain.baseName}.vcf.gz
    """
}

process SOMPY_BENCHMARK {
    // Container tag found wrong via execution 2026-08-30 (quay.io/biocontainers/hap.py:0.3.15-0
    // doesn't exist -- biocontainers tags are always "<version>--<conda-build-string>", double
    // dash, not "<version>-<N>") then fixed and confirmed with a real `docker pull`. Reused here
    // for som.py since it ships in the same package/image as hap.py -- see module header comment
    // for why som.py, not hap.py, is the right tool for a somatic (GT-less) VCF comparison.
    container 'quay.io/biocontainers/hap.py:0.3.14--py27h5c5a3ab_0'
    publishDir "${params.outdir}/benchmarking", mode: 'copy'

    input:
    tuple path(truth_vcf), path(truth_vcf_index)
    path(query_vcf)
    path(query_vcf_index)
    path(reference_fasta)
    path(reference_fai)

    output:
    path("som.stats.csv"), emit: stats
    path("som.summary.csv"), emit: summary
    path("som.metrics.json"), emit: metrics
    path("som.*"), emit: all_outputs

    script:
    // `-o som` sets the output filename prefix. Confirmed via som.py's own source (not guessed):
    // bare `-N` unconditionally writes som.stats.csv + som.metrics.json; `--happy-stats` adds the
    // friendlier som.summary.csv table read like hap.py's own summary.csv would have been.
    """
    som.py \\
        ${truth_vcf} \\
        ${query_vcf} \\
        -r ${reference_fasta} \\
        -o som \\
        -N \\
        --happy-stats
    """
}
