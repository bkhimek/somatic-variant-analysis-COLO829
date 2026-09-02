// Module 8 -- Oncogenicity/actionability interpretation layer (added 2026-09-02, after the
// real-data run -- docs/REAL_DATA_RESULTS.md -- produced 8 real PASS calls worth interpreting
// automatically rather than by hand).
//
// Data source decision: COSMIC's Cancer Gene Census (the plan's original Module 8 source) has
// had its registration on hold since Phase 0 (docs/PHASE0_FINDINGS.md) -- no account category
// cleanly fits a personal, non-commercial portfolio project. CIViC (civicdb.org) was documented
// there as the open, no-registration fallback, and is what this module actually uses.
//
// Design fork, decided with the user 2026-09-02: CIViC names variants by protein change
// ("V600E"), not genomic coordinate, so matching a Mutect2 call to CIViC evidence needs an
// intermediate protein-level annotation step. Chose full SnpEff annotation over a gene-level-only
// fallback -- more moving parts, but the only version of this module that's actually reusable on
// a future variant this project hasn't already manually looked up.
//
// Every choice below was checked by real execution, not assumed:
//   - Container tag `quay.io/biocontainers/snpeff:5.4.0c--hdfd78af_0` -- confirmed via the
//     quay.io tag API (sorted by actual last-modified date, not tag-string order) and a real
//     `docker pull`, 2026-09-02.
//   - Database choice `hg38` (SnpEff's UCSC/RefSeq-convention build), not `GRCh38.*` (SnpEff's
//     Ensembl-convention builds) -- Homo_sapiens_assembly38.fasta uses `chr`-prefixed contigs
//     (docs/data_sources.md §2), and Ensembl-convention databases use bare contig names ("1" not
//     "chr1") -- the same naming mismatch this project already hit once with the ENA GRCh37 BAMs.
//     Confirmed correct by annotating the real chr7:140753336 BRAF call and getting back
//     `NM_004333.6:c.1799T>A:p.Val600Glu` -- an exact match to Cellosaurus's documented COLO829
//     genotype and to OncoKB/ClinVar's own annotation of that exact hg38 coordinate.
//   - SnpEff's own database-hosting has a real, documented failure history (GitHub issue #596,
//     downloads failing from `snpeff.blob.core.windows.net`) -- the URLs returned by `snpEff
//     databases` in this session point at `snpeff-public.s3.amazonaws.com` instead, suggesting
//     that's been migrated, and the download succeeded cleanly in this session's real test. Not
//     re-verified beyond that; if SNPEFF_DOWNLOAD fails at runtime with a download error, that's
//     this known upstream risk resurfacing, not a bug in this pipeline.
//   - PASS-only pre-filtering before SnpEff, not an afterthought: the first real test ran SnpEff
//     directly against the full `filtered.vcf.gz` (276 candidate records, including complex
//     multiallelic sites Mutect2 filtered out for germline/contamination/strand-bias reasons) and
//     it crashed with `java.lang.OutOfMemoryError` partway through. Since this module only ever
//     needs the PASS calls, filtering first is both the fix and the correct scope -- not a
//     workaround for a limitation we're stuck with.
//   - CIViC's live GraphQL schema (`gene(entrezSymbol:)` -> `variants` -> `molecularProfiles` ->
//     `evidenceItems{evidenceLevel evidenceType significance disease{name} therapies{name}}`) was
//     confirmed by a real query against https://civicdb.org/api/graphql returning real data
//     (BRAF A598V's evidence item), not taken from CIViC's own docs, which explicitly warn
//     evidence items hang off molecular profiles rather than variants directly -- true, and the
//     query in bin/civic_annotate.py reflects it.
//   - `python:3.11-slim` (official Docker Hub image, not a `quay.io/biocontainers/*` one like the
//     rest of this pipeline) for CIVIC_ANNOTATE -- deliberate: this step is plain HTTP+JSON via
//     Python's standard library, no bioinformatics tool involved, so a bioconda-recipe container
//     would add nothing but an extra image to pull.
//
// Known, disclosed scope limit: bin/civic_annotate.py only recognizes two HGVS.p shapes --
// missense (p.Val600Glu) and frameshift (p.Ala68fs) -- because those are the only two this
// project's real Mutect2 output actually produced. A future in-frame indel or stop-gain call
// would fall through to "no protein-level match attempted" rather than a silently wrong guess;
// extending the parser is straightforward if/when a real call needs it.

process FILTER_PASS_VARIANTS {
    tag "PASS-only pre-filter"
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'   // bgzip only here -- same image already confirmed working for bgzip/tabix in modules/benchmarking.nf's PREPARE_TRUTH_VCF
    publishDir "${params.outdir}/oncogenicity", mode: 'copy'

    input:
    path(filtered_vcf)

    output:
    path("pass_only.vcf.gz"), emit: vcf

    script:
    // Plain awk on the FILTER column (field 7), not bcftools -- bcftools isn't guaranteed
    // present in every environment this pipeline runs in (confirmed missing on the user's own
    // WSL host during this module's own real-execution testing, 2026-09-02), and awk/bgzip cover
    // exactly what's needed here without adding a dependency.
    """
    zcat ${filtered_vcf} | awk '\$1 ~ /^#/ || \$7=="PASS"' | bgzip > pass_only.vcf.gz
    """
}

process SNPEFF_DOWNLOAD {
    tag "SnpEff hg38 database download"
    container 'quay.io/biocontainers/snpeff:5.4.0c--hdfd78af_0'
    publishDir "${params.snpeff_data_dir}", mode: 'copy'

    output:
    path("hg38"), emit: db_dir

    script:
    // Downloads into the process's own work dir (SnpEff's -dataDir), then publishDir copies the
    // resulting hg38/ subdirectory out to params.snpeff_data_dir -- next run's auto-detect (see
    // workflows/somatic.nf) finds it there and skips this process entirely. Confirmed working via
    // a real execution, 2026-09-02 -- ~1-2 minutes, no manual retry needed that time, though see
    // this module's header comment on SnpEff's own historical database-hosting flakiness.
    """
    snpEff download -dataDir \$(pwd) -v hg38
    """
}

process SNPEFF_ANNOTATE {
    tag "SnpEff hg38 annotation"
    container 'quay.io/biocontainers/snpeff:5.4.0c--hdfd78af_0'
    publishDir "${params.outdir}/oncogenicity", mode: 'copy'

    input:
    path(pass_vcf)
    path(snpeff_data_dir)

    output:
    path("annotated_pass.vcf"), emit: vcf

    script:
    // -dataDir must point at the directory that CONTAINS hg38/, not hg38/ itself. Nextflow
    // stages the snpeff_data_dir path input into this process's work dir under its own basename
    // ("hg38", whether it came from the pre-existing-directory branch or SNPEFF_DOWNLOAD's output
    // -- both name it that), so the work dir itself (\$(pwd)) is exactly that containing
    // directory -- matches the real command confirmed working in this session's manual testing.
    """
    snpEff -dataDir \$(pwd) hg38 ${pass_vcf} > annotated_pass.vcf
    """
}

process CIVIC_ANNOTATE {
    tag "CIViC oncogenicity annotation"
    container 'python:3.11-slim'
    publishDir "${params.outdir}/oncogenicity", mode: 'copy'

    input:
    path(annotated_vcf)

    output:
    path("civic_report.tsv"), emit: report

    script:
    """
    python3 ${projectDir}/bin/civic_annotate.py ${annotated_vcf} civic_report.tsv
    """
}
