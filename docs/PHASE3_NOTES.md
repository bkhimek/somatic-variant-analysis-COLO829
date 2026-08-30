# Phase 3 — Benchmarking (Module 5): implementation notes

**Date:** 2026-08-30
**Status:** First run found one real bug (wrong container tag), since fixed — see "First run — findings" below. Same "review isn't the same as execution" caveat as every phase so far applies.

**Scope decision made with the user before building this:** the only Mutect2 output that has actually run successfully so far is the Phase 2 `dev`-profile result — a driver-gene-restricted test against a 10,000-read-pair subsample, which came back with **zero variant calls** (documented in `PHASE2_NOTES.md`'s "Interval BED" caveat). Real full-genome Mutect2 execution is still deferred (needs interval-scatter/gather architecture plus likely another cloud compute burst — `PHASE2_NOTES.md`'s sign-off section). Rather than block Module 5 on that unfinished work, the decision was to **build and wire in the benchmarking module now, and validate it against the existing (empty) `dev`-profile output purely to prove the DAG/container/hap.py command line are correct** — exactly the same "code correctness before real signal" split that Phase 1 (chr21 test) and Phase 2 (driver-gene BED test) already established. Real precision/recall/F1 numbers stay deferred until real full-genome Mutect2 output exists.

**What that means concretely:** running this against current data should produce a `happy.summary.csv` showing 0 true positives, 0 false positives, and however many false negatives the truth set contains (since there's essentially nothing in the query VCF to match). That's the expected, correct outcome for this test — not a bug. Do not read anything into the actual precision/recall numbers from this run; the only thing worth checking is that hap.py runs to completion without erroring.

---

## First run — findings (2026-08-30)

**Real bug: `PREPARE_TRUTH_VCF`'s container tag, `quay.io/biocontainers/htslib:1.21--h566b1c6_0`, doesn't exist** — Docker failed with `manifest unknown`. That tag was pattern-matched off the already-confirmed `samtools:1.21--h50ea8bc_0` image's build-hash convention rather than actually checked against the registry, which was the mistake — this project's own established standard (see the bwa-mem2 memory correction, the Mutect2 memory research) is to verify container/version claims for real before shipping them, and this one slipped through. Verified properly via research afterwards rather than guessing a second time: real htslib 1.21 tags do exist (`1.21--h566b1c6_1`, note the `_1` not `_0`, or `1.21--h5efdd21_0`), but the better fix is to **not pull a second image at all** — `samtools:1.21--h50ea8bc_0` (already confirmed working in this pipeline since Phase 1) bundles `bgzip`/`tabix` too, because bioconda's samtools recipe pulls in htslib as a genuine runtime dependency (via htslib's own `run_exports` mechanism, not just a build-time link), which installs those binaries into the same conda environment. **Fix:** `PREPARE_TRUTH_VCF` now reuses `samtools:1.21--h50ea8bc_0` instead of a separate htslib image.

**Cosmetic, non-blocking: a Nextflow `WARN` about a useless `.first()` call.** `PREPARE_TRUTH_VCF.out.truth_vcf_indexed.first()` triggered `WARN: The operator 'first' is useless when applied to a value channel which returns a single value by definition`. Harmless (the channel already only ever emits once either way) but easy to clean up: `PREPARE_TRUTH_VCF`'s only input is a plain value (not a channel), so Nextflow already treats a process with all-singleton-value inputs as producing a value-channel output — `.first()` on it is redundant. Removed. (The identical latent pattern exists in `INDEX_FASTA`/`CREATE_SEQUENCE_DICTIONARY`'s `.first()` calls from Phase 2, but their build branches have never actually executed in any run so far — always hit the pre-existing-`.fai`/`.dict` shortcut instead — so this hasn't surfaced there yet. Left alone rather than touched, since that's already-signed-off Phase 2 code and this warning is cosmetic, not a correctness issue.)

---

## What's implemented

- `modules/benchmarking.nf` — `PREPARE_TRUTH_VCF` (bgzip + tabix on NYGC's plain uncompressed truth VCF), `HAPPY_BENCHMARK` (hap.py comparing the truth VCF against `FILTER_MUTECT_CALLS`' `filtered.vcf.gz`)
- `workflows/somatic.nf` — extended (not forked) with Module 5, per the established "each phase extends the same workflow" design; new `truth_set_vcf` take, auto-detect-or-build for the bgzip+tabix truth VCF (same pattern as the `.fai`/`.dict` auto-prep from Phase 2), new `happy_summary` emit
- `main.nf` — `truth_set_vcf` added to the required-params fail-fast check; passed through to `SOMATIC`
- `nextflow.config` — `truth_set_vcf`/`truth_set_bed` params already existed as Phase 0 scaffolding (`null` placeholders) — no structural change needed, just updated the comment now that they're actually wired in and confirmed no BED restriction applies

---

## Design choices made, and why (so you can sanity-check them)

1. **hap.py's default comparison engine (xcmp), not vcfeval.** vcfeval (Real Time Genomics' engine) is generally considered more accurate for complex indels, but needs its own SDF-format reference build and the `rtg-tools` package as an extra dependency. xcmp is hap.py's own default, needs nothing beyond what's already in the `hap.py` container, and is a perfectly standard choice for a first benchmarking pass. Documented here as a real choice, not an oversight — worth revisiting once real full-genome results exist and indel accuracy specifically becomes worth scrutinising closely.
2. **No confidence-region BED restriction.** Confirmed in Phase 0 (`data_sources.md` §3): NYGC does not publish a high-confidence callable-regions BED for SNV/indel, only a CNV bed and SV bedpe (neither applicable here). hap.py runs unrestricted across the whole genome. This was a Phase 0 finding, not a new decision.
3. **No `--pass-only` flag needed.** hap.py's own default already only counts `PASS` variants as positive calls on both the truth and query side. `FILTER_MUTECT_CALLS`' `filtered.vcf.gz` marks non-PASS calls with real filter reasons (contamination, orientation bias, etc.) — hap.py's default behaviour of excluding those from being counted as calls is exactly what we want, not something to override.
4. **Truth VCF bgzip/tabix handled automatically, mirroring Phase 2's resource-VCF pattern**, rather than asking you to prep it by hand first. `data_sources.md` already confirmed the file NYGC ships is plain uncompressed VCF, not `.vcf.gz` — `PREPARE_TRUTH_VCF` bgzips and indexes it the first time, and auto-detects/skips on later runs if you've already done this yourself.

---

## Real, unverified risk worth flagging before you run this

**Contig-naming and reference-build match between the truth VCF and our reference FASTA has not been directly checked yet.** The NYGC truth set was generated by NYGC's own pipeline against their own reference build (`data_sources.md` §3's lineage-mismatch caveat already covers this at the calling-methodology level) — but the more basic, mechanical question of whether its VCF header uses `chr`-prefixed contig names (matching `Homo_sapiens_assembly38.fasta`, the Broad reference this pipeline aligns against, per `dev_intervals.bed`'s note that Broad's convention is `chr`-prefixed) hasn't been confirmed by actually looking at the truth VCF's header. If they don't match, hap.py will fail fast with a clear contig-mismatch error rather than silently producing wrong numbers, so this isn't a silent-corruption risk — but it's worth checking (`zcat truth_set/COLO-829--COLO-829BL.snv.indel.final.v6.annotated.vcf.gz | grep '^##contig' | head` once bgzipped, or `head -50` on the plain VCF before that) if hap.py errors out on the first real run, before assuming it's something else.

---

## Things you need to check/fix before this actually runs

1. **This has never been executed.** Same expectation as every phase so far.
2. **`hap.py` container tag was verified to exist in Phase 0 but never actually pulled/run** (`quay.io/biocontainers/hap.py:0.3.15-0`) — confirm `docker pull` works on your machine before the full pipeline run, same as the samtools/GATK images were spot-checked in earlier phases.
3. ~~`htslib` container tag for `PREPARE_TRUTH_VCF` hasn't been confirmed yet~~ — **found wrong via execution, fixed:** see "First run — findings" above. `PREPARE_TRUTH_VCF` now reuses `samtools:1.21--h50ea8bc_0`, already confirmed working since Phase 1, so no new image needs pulling for this step.
4. **Memory hasn't been preemptively tuned for either new process** — both inherit the flat profile default. Given this project has now hit two real GATK-in-Docker/JVM memory surprises in Phase 2, this is flagged rather than assumed fine: hap.py on a genome-wide (if variant-free) comparison against a ~36MB truth VCF is expected to be lightweight, but if it OOMs, that's a real finding to fix with evidence (same standard applied throughout this project), not a reason to have guessed a bigger number upfront.
5. **`docs/benchmarking_results.md`** (referenced in the repo layout, not created yet) stays unwritten until there's a real result worth recording — writing it now with a 0/0/0 placeholder table would just need to be redone once real full-genome numbers exist. When that happens, carry forward the two documented caveats: the NYGC lineage-mismatch note (`data_sources.md` §3 — frame results as "agreement with a GATK-family-adjacent somatic pipeline, itself independently validated against Craig et al. 2016 at 98% concordance" rather than unmediated ground truth) and the no-confidence-region-BED caveat above.

---

## How to run it (once you've checked the items above)

```bash
conda activate nextflow
cd ~/projects/somatic-variant-analysis-COLO829

nextflow run main.nf -profile docker,dev \
    --tumour_reads_1        fastq_dev_sample/TUMOUR_R1.sample.fastq.gz \
    --tumour_reads_2        fastq_dev_sample/TUMOUR_R2.sample.fastq.gz \
    --normal_reads_1        fastq_dev_sample/NORMAL_R1.sample.fastq.gz \
    --normal_reads_2        fastq_dev_sample/NORMAL_R2.sample.fastq.gz \
    --reference_fasta       reference/Homo_sapiens_assembly38.fasta \
    --panel_of_normals      reference/1000g_pon.hg38.vcf.gz \
    --germline_resource     reference/af-only-gnomad.hg38.vcf.gz \
    --common_biallelic_sites reference/small_exac_common_3.hg38.vcf.gz \
    --truth_set_vcf         truth_set/COLO-829--COLO-829BL.snv.indel.final.v6.annotated.vcf \
    -resume
```

(Path above for `--truth_set_vcf` matches where `data_sources.md` §3 says the file landed — adjust if you put it elsewhere.) `-resume` reuses everything already computed for Modules 1–4, so this should only actually execute the new `PREPARE_TRUTH_VCF` + `HAPPY_BENCHMARK` steps.

Expected new outputs:
- `results/benchmarking/truth_set/<truth_vcf_basename>.vcf.gz` (+ `.tbi`)
- `results/benchmarking/happy.summary.csv` — the headline precision/recall/F1 table (expect all zeros/near-zeros against current dev-profile data — see "Status" above)
- `results/benchmarking/happy.extended.csv`, plus hap.py's other standard output files (`happy.vcf.gz`, ROC curve data, etc.)

Per the caveat above, don't be surprised (or alarmed) by a summary showing 0 TP/0 FP against the current data — check that the pipeline **ran to completion without errors** first, same split every phase so far has established.
