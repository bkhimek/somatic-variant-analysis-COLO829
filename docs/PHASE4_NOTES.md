# Phase 4 — Copy-number calling (Module 6): implementation notes

**Date:** 2026-08-31
**Status:** Built, not yet run. Same "review isn't the same as execution" caveat as every phase so far applies — Phases 1-3 each found real bugs only by actually executing the code, and there's no reason to expect this phase is different.

**How this phase got scoped:** unlike Phases 1-3, the project plan's exact CNVkit section (module number, specific flags/resources it names) wasn't retrieved before building this — the user explicitly asked for it to be designed from CNVkit's standard tumour/normal WGS workflow instead of blocking on that (2026-08-31). Everything below was researched against CNVkit's own current documentation and a couple of independently-verified facts (see "Design choices" below), not guessed — but if the plan specifies something different (a particular module number, a required annotation file, specific flags), that should override what's here.

---

## What's implemented

- `modules/cnvkit.nf` — `CNVKIT_BATCH` (CNVkit's whole-genome (`-m wgs`) tumour/normal workflow: builds its own reference from the normal BAM, computes bins against the reference FASTA's accessible regions, produces `.cnr`/`.cns`), `CNVKIT_CALL` (converts segments to integer copy-number calls)
- `workflows/somatic.nf` — extended (not forked) with Module 6, reusing the same `tumour_bam_ch`/`normal_bam_ch` tuples Module 4 already branches out of `MARK_DUPLICATES.out.bam` — no new branching logic needed; new `cnvkit_cnr`/`cnvkit_call` emits
- `main.nf` — no new required params (CNVkit only needs `reference_fasta`, already required since Phase 1, and the dedup BAMs, already produced by Module 2)
- `nextflow.config` — no changes; CNVkit inherits the flat per-process default (4 cpus / 8GB) like most processes in this pipeline, see "Things to check" below for why that's flagged rather than assumed fine

---

## Design choices made, and why (so you can sanity-check them)

1. **`-m wgs` (whole-genome mode), no target/antitarget BED.** Confirmed via CNVkit's own docs (`nonhybrid.html`): in WGS mode, the reference genome's sequencing-accessible regions ("access" BED) are used as CNVkit's "targets" and are computed on the fly if not supplied — there's no separate access-BED prep step needed, unlike a capture/exome workflow.
2. **No auto-detect-or-build wrapper, unlike `.fai`/`.dict`/the truth VCF.** `INDEX_FASTA`, `CREATE_SEQUENCE_DICTIONARY`, and `PREPARE_TRUTH_VCF` all follow an "auto-detect existing file, else build it" pattern because those artifacts are expensive-ish, shared, and reused identically across every run. CNVkit's own internal reference-building (from the normal BAM) and access-BED computation are comparatively cheap, per-run, and tied to whatever tumour/normal pair is being analysed that run — there's nothing worth caching outside the process, so `CNVKIT_BATCH` just does it inline every time, matching CNVkit's own designed usage pattern rather than working against it.
3. **`-y`/`--male-reference` included, verified rather than guessed.** Getting sex-chromosome ploidy handling wrong silently biases every chrX/chrY copy-number call (CNVkit defaults to assuming a female/XX reference otherwise). Checked before writing any code, not assumed: [ATCC's own COLO 829 (CRL-1974) product page](https://www.atcc.org/products/crl-1974) lists the donor as a 45-year-old male. `-y` is passed to both `CNVKIT_BATCH` (which builds the reference internally) and `CNVKIT_CALL` — CNVkit's own docs are explicit that the same flag must be repeated at `call` time for consistent handling, not just at reference-build time.
4. **`--drop-low-coverage` included.** CNVkit's own documentation recommends this specifically "for tumor samples" — it drops bins at/near zero read depth from segmentation rather than treating them as real (very negative) log2-ratio signal. This matters more than usual here: the only tumour BAM that exists so far is Phase 2's `dev`-profile result, a 10,000-read-pair subsample spread across the whole genome, so the overwhelming majority of whole-genome bins will have exactly zero coverage.
5. **No `--annotate refFlat.txt`, deferred.** CNVkit's docs describe a gene-annotation database as needed for *readable per-target gene labels* in the output, not for `batch` to run at all. Sourcing, downloading, and verifying a GRCh38 refFlat.txt is a real task (another resource to pin a version/URL for in `docs/data_sources.md`, in the same spirit as the GATK resource bundles) — deferred rather than rushed in for a first cut whose query BAM has almost no real coverage to label meaningfully anyway. Outputs will carry generic/positional bin labels instead of gene symbols until this is revisited.
6. **No `--scatter`/`--diagram`, deferred.** Both are valid, well-documented `batch` flags for QC plots. Deliberately not added this time — Phase 3 just got burned by adding `som.py --happy-stats` (a flag that turned out to need companion flags this project wasn't ready to pick correctly) without checking what it actually required to succeed first. Same risk shape here: better to prove the base `.cnr`/`.cns`/`.call.cns` path runs cleanly first, then add plotting flags once that's confirmed, rather than risk a repeat of exactly the mistake just made and documented in `docs/PHASE3_NOTES.md`.
7. **Two processes (`CNVKIT_BATCH` then `CNVKIT_CALL`), not one.** Confirmed via CNVkit's own pipeline documentation: `batch` does NOT run the `call` step automatically — it stops at segmented-but-not-yet-integer `.cns` output. `call` is a genuinely separate step, so it gets its own process, consistent with this pipeline's one-tool-invocation-per-process style elsewhere (e.g. Mutect2 → LearnReadOrientationModel → FilterMutectCalls as three processes, not one).

---

## Real, unverified risks worth flagging before you run this

1. **Container tag not independently verified.** `quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0` is the plan's original choice, and Phase 0 confirmed the underlying bioconda recipe (cnvkit 0.9.10, build 0) exists — but the exact quay.io tag *string* was never pulled and confirmed, because quay.io blocks automated tag-list fetching from this sandbox (robots.txt) the same way it did for htslib and hap.py, both of which turned out to have wrong guessed tags (see `docs/PHASE3_NOTES.md`'s "First run" and "Second run" findings). Given that exact history repeating twice already, **please run `docker pull quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0` yourself and confirm it succeeds before wiring this into a real run** — same process that worked for confirming the corrected hap.py tag in Phase 3. If it fails, don't guess a replacement blind; tell me the exact error and I'll research a verified alternative the same way.
2. **Whether CNVkit tolerates an almost-entirely-zero-coverage genome-wide BAM is unverified.** `--drop-low-coverage` is the documented mitigation, but this project's `dev`-profile tumour BAM is a more extreme case than typical "some bins are low-coverage" tumour heterogeneity — it's a 10,000-read-pair subsample against the *whole genome*, meaning close to zero real coverage almost everywhere. It's plausible CNVkit's `fix`/`segment` steps handle "nearly all bins dropped" gracefully (similar to how Mutect2 and som.py both handled "structurally nothing to find" gracefully in Phases 2-3) — but it's also plausible something in that chain (division by a near-zero mean, a segmentation algorithm given almost no retained bins) errors out instead. This is a real, unresolved question this project's usual "find out by executing" approach will answer, not something assumed away.
3. **Memory/CPU not preemptively tuned.** CNVkit inherits the flat per-process default (4 cpus / 8GB from the active profile). Given this project has now hit real GATK-in-Docker and bwa-mem2 memory surprises more than once (Phases 1-2), this is flagged rather than assumed fine — a first `dev`-profile run against a near-empty BAM is expected to be lightweight regardless of how it behaves, so this risk mostly matters once a real `full`-profile, real-coverage BAM exists.
4. **Runtime of the on-the-fly access-BED computation across the whole GRCh38 reference is unmeasured.** CNVkit computes this by scanning the reference FASTA for N-gaps every time `batch` runs (no caching, per design choice #2 above) — typically a few minutes for a human genome per public reports, but not independently timed in this project yet.

---

## How to run it (once you've checked the items above)

```bash
conda activate nextflow
cd ~/projects/somatic-variant-analysis-COLO829

# Confirm the container tag first (see "Real, unverified risks" #1 above) --
# if this fails, stop and report back rather than guessing a replacement tag.
docker pull quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0

nextflow run main.nf -profile docker,dev \
    --tumour_reads_1        fastq_dev_sample/TUMOUR_R1.sample.fastq.gz \
    --tumour_reads_2        fastq_dev_sample/TUMOUR_R2.sample.fastq.gz \
    --normal_reads_1        fastq_dev_sample/NORMAL_R1.sample.fastq.gz \
    --normal_reads_2        fastq_dev_sample/NORMAL_R2.sample.fastq.gz \
    --reference_fasta       reference/Homo_sapiens_assembly38.fasta \
    --panel_of_normals      reference/1000g_pon.hg38.vcf.gz \
    --germline_resource     reference/af-only-gnomad.hg38.vcf.gz \
    --common_biallelic_sites reference/small_exac_common_3.hg38.vcf.gz \
    --truth_set_vcf         truth_set/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf \
    -resume
```

`-resume` reuses everything already computed through Modules 1-5 (Phase 3 left that fully cached and signed off) — this should only actually execute the new `CNVKIT_BATCH` and `CNVKIT_CALL` steps.

Expected new outputs:
- `results/cnvkit/*.cnr` — bin-level log2 copy ratios
- `results/cnvkit/*.cns` — segmented (but not yet integer) copy-number regions
- `results/cnvkit/*.call.cns` — integer copy-number calls (the clinically-meaningful output)
- `results/cnvkit/*.cnn` — coverage files CNVkit uses internally to build its per-run reference

**What a successful run here would prove:** the container pulls and runs, CNVkit's WGS tumour/normal batch pipeline executes against real BAMs and a real reference without erroring, and `call` correctly consumes `batch`'s segment output. **What it would NOT prove:** anything about COLO829's actual copy-number profile — same "code correctness, not real signal" split established in every phase so far, since the only tumour BAM that exists is the near-empty `dev`-profile subsample. Expect `.call.cns` to show mostly (or entirely) neutral/diploid calls, not real amplifications or deletions.

Don't be surprised (or alarmed) by a flat, uninteresting `.call.cns` — check that the pipeline **ran to completion without errors** first, same split every phase so far has established. If `CNVKIT_BATCH` or `CNVKIT_CALL` fails outright, that's the more interesting and more likely finding (see "Real, unverified risks" #2 above) — paste the error back and we'll diagnose it the same way every other real bug in this project has been diagnosed: from the actual failure, not a guess.
