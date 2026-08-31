// Phase 1: QC (Module 1) + Alignment (Module 2). Phase 2 (added 2026-08-30): Contamination
// estimation (Module 3) + Mutect2 somatic calling (Module 4). Tumour and normal are processed
// independently through Modules 1-2 and Module 3's pileup step (per plan §6 Module 2 note --
// "Same module runs twice"), then rejoined for CalculateContamination and Mutect2 itself, which
// genuinely need both samples together. Phase 3 (added 2026-08-30): Benchmarking (Module 5) --
// compares FILTER_MUTECT_CALLS' output against the NYGC COLO829 truth set via som.py (NOT
// hap.py -- see modules/benchmarking.nf's header comment for why a GT-less somatic VCF needs
// som.py specifically). See docs/PHASE3_NOTES.md for what this run can and can't validate yet
// (short version: the only Mutect2 output that exists so far has zero variants, so this proves
// the DAG/container/som.py command line work, not real precision/recall numbers).
//
// Phase 4 (added 2026-08-31): Copy-number calling (Module 6) -- CNVkit's whole-genome tumour/
// normal workflow against the same dedup BAMs Module 4 already uses, then integer copy-number
// calls. See modules/cnvkit.nf's header comment for scope limits (no --annotate/--scatter/
// --diagram yet) and docs/PHASE4_NOTES.md for what this run can and can't validate.
//
// Phase 5+ (SigProfiler, interpretation, ...) will extend this same workflow file rather than
// starting a new one, so the DAG stays in one place end to end.

include { FASTQC; MULTIQC }                                             from '../modules/fastqc.nf'
include { BWA_MEM2_INDEX; BWA_MEM2_ALIGN; SAMTOOLS_SORT; MARK_DUPLICATES } from '../modules/alignment.nf'
include { INDEX_FASTA; CREATE_SEQUENCE_DICTIONARY; INDEX_VCF }           from '../modules/reference_prep.nf'
include { GET_PILEUP_SUMMARIES; CALCULATE_CONTAMINATION }                from '../modules/contamination.nf'
include { MUTECT2; LEARN_READ_ORIENTATION_MODEL; FILTER_MUTECT_CALLS }   from '../modules/mutect2.nf'
include { PREPARE_TRUTH_VCF; SOMPY_BENCHMARK }                          from '../modules/benchmarking.nf'
include { CNVKIT_BATCH; CNVKIT_CALL }                                   from '../modules/cnvkit.nf'

workflow SOMATIC {

    take:
    samples_ch              // channel: tuple(sample_id, reads1, reads2) — one entry per sample (tumour, normal)
    reference_fasta          // path: reference genome FASTA
    panel_of_normals         // path: Mutect2 --panel-of-normals VCF (gz)
    germline_resource        // path: Mutect2 --germline-resource VCF (gz)
    common_biallelic_sites   // path: GetPileupSummaries -V VCF (gz)
    truth_set_vcf            // path: NYGC COLO829 truth-set VCF (plain, uncompressed)

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
    // building it, but NOT the ~19GB an earlier version of this comment claimed -- that figure
    // was an unverified guess and turned out wrong; see docs/PHASE1_NOTES.md's 2026-08-30 update
    // for real-world numbers (commonly 30GB+, per bwa-mem2's own GitHub issues) and the actual
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

    // ---- Reference/resource prep for Modules 3-4 (added 2026-08-30) ----
    // .fai/.dict: flagged as an open TODO since Phase 0 (docs/PHASE0_FINDINGS.md §7 action 8),
    // never actually built until now. Auto-detect/build, same pattern as the bwa-mem2 index
    // above -- these are cheap (seconds), not a repeat of that memory saga, but still worth
    // skipping if you already pulled them from the Broad FTP per docs/data_sources.md §2.
    def fai_marker = file("${reference_fasta}.fai")
    if (fai_marker.exists()) {
        log.info "Existing .fai found alongside ${reference_fasta} -- skipping INDEX_FASTA"
        fai_ch = Channel.value(fai_marker)
    } else {
        INDEX_FASTA(reference_fasta)
        fai_ch = INDEX_FASTA.out.fai.first()
    }

    def dict_marker = file("${reference_fasta}".replaceAll(/\.(fa|fasta|fna)$/, '') + '.dict')
    if (dict_marker.exists()) {
        log.info "Existing .dict found alongside ${reference_fasta} -- skipping CREATE_SEQUENCE_DICTIONARY"
        dict_ch = Channel.value(dict_marker)
    } else {
        CREATE_SEQUENCE_DICTIONARY(reference_fasta)
        dict_ch = CREATE_SEQUENCE_DICTIONARY.out.dict.first()
    }

    // Resource VCF indices (.tbi). Unlike the bwa-mem2 index and .fai/.dict above, this does
    // NOT auto-detect a pre-existing .tbi and skip -- deliberately simpler, for two reasons:
    // (1) `gatk IndexFeatureFile` takes seconds per VCF, so there's no expensive-rebuild problem
    // like bwa-mem2's to optimise away; (2) Nextflow only cleanly supports invoking the same
    // process multiple times in one workflow if each call's return value is captured into its
    // own variable rather than read back via the shared, ambiguous `PROCESS.out` after more than
    // one call -- simplest and safest is one call over a combined channel of all three VCFs,
    // then split the results back apart below by filename (all three have distinct basenames).
    // -resume still skips this entirely on every run after the first regardless.
    resource_vcfs_ch = Channel.fromPath([panel_of_normals, germline_resource, common_biallelic_sites])
    INDEX_VCF(resource_vcfs_ch)

    indexed_resources_branched = INDEX_VCF.out.vcf_with_index.branch {
        pon:      it[0].getName() == file(panel_of_normals).getName()
        germline: it[0].getName() == file(germline_resource).getName()
        common:   it[0].getName() == file(common_biallelic_sites).getName()
    }
    panel_of_normals_ch       = indexed_resources_branched.pon.first()
    germline_resource_ch      = indexed_resources_branched.germline.first()
    common_biallelic_sites_ch = indexed_resources_branched.common.first()

    // NO_FILE sentinel (see assets/NO_FILE) -- params.interval_list is null on the `full`
    // profile (no restriction); on `dev` it's data/gene_lists/dev_intervals.bed
    // (docs/PHASE2_NOTES.md has the important caveat on what that restriction does and
    // doesn't validate against the current dev-profile subsample BAMs).
    interval_ch = params.interval_list
        ? Channel.value(file(params.interval_list))
        : Channel.value(file("${projectDir}/assets/NO_FILE"))

    // ---- Module 3: Contamination estimation ----
    // GET_PILEUP_SUMMARIES runs on both samples generically (like Module 1/2) -- .branch{}
    // below is what actually knows which sample is tumour vs normal, needed because
    // CalculateContamination and Mutect2 (unlike Modules 1-2) require tumour/normal paired
    // together, not processed identically-and-independently.
    GET_PILEUP_SUMMARIES(
        MARK_DUPLICATES.out.bam,
        reference_fasta, fai_ch, dict_ch,
        common_biallelic_sites_ch,
        interval_ch
    )

    dedup_bams_branched = MARK_DUPLICATES.out.bam.branch {
        tumour: it[0] == 'COLO829'
        normal: it[0] == 'COLO829BL'
    }
    tumour_bam_ch = dedup_bams_branched.tumour
    normal_bam_ch = dedup_bams_branched.normal

    pileups_branched = GET_PILEUP_SUMMARIES.out.pileups.branch {
        tumour: it[0] == 'COLO829'
        normal: it[0] == 'COLO829BL'
    }

    CALCULATE_CONTAMINATION(
        pileups_branched.tumour.map { sample_id, table -> table },
        pileups_branched.normal.map { sample_id, table -> table }
    )

    // ---- Module 4: Mutect2 somatic calling ----
    MUTECT2(
        tumour_bam_ch, normal_bam_ch,
        reference_fasta, fai_ch, dict_ch,
        panel_of_normals_ch, germline_resource_ch,
        interval_ch
    )

    LEARN_READ_ORIENTATION_MODEL(MUTECT2.out.f1r2)

    FILTER_MUTECT_CALLS(
        MUTECT2.out.vcf, MUTECT2.out.vcf_index, MUTECT2.out.stats,
        reference_fasta, fai_ch, dict_ch,
        CALCULATE_CONTAMINATION.out.contamination_table,
        CALCULATE_CONTAMINATION.out.segments_table,
        LEARN_READ_ORIENTATION_MODEL.out.model
    )

    // ---- Module 5: Benchmarking (Phase 3, added 2026-08-30) ----
    // Auto-detect-or-build for the truth VCF's bgzip+tabix companion, same pattern as .fai/.dict
    // above -- skip PREPARE_TRUTH_VCF if you've already bgzipped/indexed it yourself alongside
    // the plain VCF NYGC ships.
    def truth_vcf_gz_marker  = file("${truth_set_vcf}.gz")
    def truth_vcf_tbi_marker = file("${truth_set_vcf}.gz.tbi")
    if (truth_vcf_gz_marker.exists() && truth_vcf_tbi_marker.exists()) {
        log.info "Pre-built bgzip+tabix truth VCF found alongside ${truth_set_vcf} -- skipping PREPARE_TRUTH_VCF"
        truth_vcf_ch = Channel.value([truth_vcf_gz_marker, truth_vcf_tbi_marker])
    } else {
        // No .first() here (unlike INDEX_FASTA/CREATE_SEQUENCE_DICTIONARY above) -- found via
        // execution 2026-08-30: PREPARE_TRUTH_VCF's only input (truth_set_vcf) is a plain value,
        // not a channel, so Nextflow already treats its output as a value channel (a process
        // whose every input is a singleton value runs once and emits a value channel by
        // definition). Calling .first() on it is harmless but triggers a "useless" WARN --
        // removed here since it's simple to avoid. The same latent (harmless) warning would
        // apply to INDEX_FASTA/CREATE_SEQUENCE_DICTIONARY's .first() calls above if their
        // build branches ever actually ran instead of hitting the pre-existing-file shortcut,
        // but touching that already-signed-off Phase 2 code isn't in scope here.
        PREPARE_TRUTH_VCF(truth_set_vcf)
        truth_vcf_ch = PREPARE_TRUTH_VCF.out.truth_vcf_indexed
    }

    SOMPY_BENCHMARK(
        truth_vcf_ch,
        FILTER_MUTECT_CALLS.out.vcf, FILTER_MUTECT_CALLS.out.vcf_index,
        reference_fasta, fai_ch
    )

    // ---- Module 6: Copy-number calling (Phase 4, added 2026-08-31) ----
    // Reuses the same tumour_bam_ch/normal_bam_ch tuples Module 4 already branched out above --
    // no new branching needed, CNVkit just wants the same dedup BAMs Mutect2 does.
    CNVKIT_BATCH(
        tumour_bam_ch, normal_bam_ch,
        reference_fasta, fai_ch,
        params.cnvkit_target_avg_size
    )

    CNVKIT_CALL(CNVKIT_BATCH.out.cns)

    emit:
    dedup_bams        = MARK_DUPLICATES.out.bam       // tuple(sample_id, bam, bai)
    dedup_metrics     = MARK_DUPLICATES.out.metrics
    multiqc_report    = MULTIQC.out.report
    contamination_table = CALCULATE_CONTAMINATION.out.contamination_table
    filtered_vcf      = FILTER_MUTECT_CALLS.out.vcf    // Phase 3's benchmarking-against-truth-set target
    sompy_stats       = SOMPY_BENCHMARK.out.stats      // TP/FP/FN counts by variant type (som.stats.csv -- no --happy-stats summary table, see modules/benchmarking.nf)
    cnvkit_cnr        = CNVKIT_BATCH.out.cnr           // bin-level log2 copy-ratio (Phase 4)
    cnvkit_call       = CNVKIT_CALL.out.call_cns       // integer copy-number segment calls (Phase 4)
}
