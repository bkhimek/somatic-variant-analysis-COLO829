// Phase 1 workflow: QC (Module 1) + Alignment (Module 2), tumour and normal each
// processed independently through the identical chain, per plan §6 Module 2 note
// ("Same module runs twice (tumour, normal) — reuse from Project 4 design").
//
// Phase 2+ (contamination, Mutect2, benchmarking, ...) will extend this same workflow
// file rather than starting a new one, so the DAG stays in one place end to end.

include { FASTQC; MULTIQC }                              from '../modules/fastqc.nf'
include { BWA_MEM2_INDEX; BWA_MEM2_ALIGN; SAMTOOLS_SORT; MARK_DUPLICATES } from '../modules/alignment.nf'

workflow SOMATIC {

    take:
    samples_ch          // channel: tuple(sample_id, reads1, reads2) — one entry per sample (tumour, normal)
    reference_fasta      // path: reference genome FASTA

    main:
    // ---- Module 1: QC ----
    FASTQC(samples_ch)

    // MultiQC aggregates across BOTH samples in one report (plan's explicit Phase 1 addition
    // beyond Project 4's single-sample QC) -- collect every FastQC zip/html before running it.
    // NOTE: FASTQC's html/zip outputs are each a List (glob-matched, 2 files per sample --
    // one per read of the pair), so this must concatenate those lists (html + zip), not nest
    // them ([html, zip]) -- nesting would hand MULTIQC a list-of-lists instead of flat files.
    fastqc_files_ch = FASTQC.out.reports
        .flatMap { sample_id, html, zip -> html + zip }
        .collect()

    MULTIQC(fastqc_files_ch)

    // ---- Module 2: Alignment ----
    // Pre-built-index auto-detection (added 2026-08-30, after the full-genome bwa-mem2 index
    // was built on a one-off high-RAM AWS instance -- see docs/data_sources.md §2 and
    // docs/PHASE1_NOTES.md for why: bwa-mem2 index needs ~87-110GB RAM to build, which no
    // desktop/laptop in this project has. Using a pre-built index for alignment needs less than
    // building it, but NOT the ~19GB this comment originally claimed -- that figure was an
    // unverified guess and turned out wrong; see docs/PHASE1_NOTES.md's 2026-08-30 update for
    // real-world numbers (commonly 30GB+, per bwa-mem2's own GitHub issues) and the actual
    // memory/cpus settings this required in nextflow.config). If the five
    // bwa-mem2 index files already sit alongside reference_fasta (same convention bwa-mem2
    // itself uses -- sibling files named "<fasta>.0123", "<fasta>.amb", etc.), skip rebuilding
    // entirely and feed BWA_MEM2_ALIGN the existing files. Otherwise (e.g. the chr21 dev-profile
    // slice, which has no pre-built index) fall back to building it in-pipeline as before --
    // Nextflow's -resume still caches that build across runs either way.
    def prebuilt_index_marker = file("${reference_fasta}.bwt.2bit.64")
    if (prebuilt_index_marker.exists()) {
        log.info "Pre-built bwa-mem2 index found alongside ${reference_fasta} -- skipping BWA_MEM2_INDEX"
        index_files_ch = Channel.fromPath("${reference_fasta}.{0123,amb,ann,bwt.2bit.64,pac}").collect()
    } else {
        // Index is built once from the reference and reused for both samples -- Nextflow
        // caches this on -resume, so it only actually runs on the first invocation.
        BWA_MEM2_INDEX(reference_fasta)

        // .first() turns the index output into a broadcastable value: BWA_MEM2_INDEX runs once
        // and emits once, but samples_ch has two items (tumour, normal). Without .first(), Nextflow
        // zips the two channels positionally and the index channel runs dry after the first sample --
        // BWA_MEM2_ALIGN would silently only fire once instead of twice.
        index_files_ch = BWA_MEM2_INDEX.out.index.first()
    }

    BWA_MEM2_ALIGN(samples_ch, reference_fasta, index_files_ch)
    SAMTOOLS_SORT(BWA_MEM2_ALIGN.out.sam)
    MARK_DUPLICATES(SAMTOOLS_SORT.out.bam)

    emit:
    dedup_bams   = MARK_DUPLICATES.out.bam       // tuple(sample_id, bam, bai) -- Phase 2 (contamination, Mutect2) consumes this
    dedup_metrics = MARK_DUPLICATES.out.metrics
    multiqc_report = MULTIQC.out.report
}
