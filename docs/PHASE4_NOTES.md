# Phase 4 — Copy-number calling (Module 6): implementation notes

**Date:** 2026-08-31
**Status:** First run found a real bug — not the container tag (that pulled fine), but a genuine statistical floor: `CNVKIT_BATCH` exited 0 without ever writing a `.cns` file, because the dev-profile subsample is too sparse for CNVkit's whole-genome binning approach to compute valid per-bin statistics anywhere. Root-caused via CNVkit's own `fix.py` source (not guessed) and fixed with a `dev`-profile-only bin-size override — see "First run — findings" below. Not yet re-run against that fix. Same "review isn't the same as execution" caveat as every phase so far applies.

**How this phase got scoped:** unlike Phases 1-3, the project plan's exact CNVkit section (module number, specific flags/resources it names) wasn't retrieved before building this — the user explicitly asked for it to be designed from CNVkit's standard tumour/normal WGS workflow instead of blocking on that (2026-08-31). Everything below was researched against CNVkit's own current documentation and a couple of independently-verified facts (see "Design choices" below), not guessed — but if the plan specifies something different (a particular module number, a required annotation file, specific flags), that should override what's here.

---

## What's implemented

- `modules/cnvkit.nf` — `CNVKIT_BATCH` (CNVkit's whole-genome (`-m wgs`) tumour/normal workflow: builds its own reference from the normal BAM, computes bins against the reference FASTA's accessible regions, produces `.cnr`/`.cns`), `CNVKIT_CALL` (converts segments to integer copy-number calls)
- `workflows/somatic.nf` — extended (not forked) with Module 6, reusing the same `tumour_bam_ch`/`normal_bam_ch` tuples Module 4 already branches out of `MARK_DUPLICATES.out.bam` — no new branching logic needed; new `cnvkit_cnr`/`cnvkit_call` emits
- `main.nf` — no new required params (CNVkit only needs `reference_fasta`, already required since Phase 1, and the dedup BAMs, already produced by Module 2)
- `nextflow.config` — CNVkit inherits the flat per-process default (4 cpus / 8GB) like most processes in this pipeline, see "Real, unverified risks" below for why that's flagged rather than assumed fine. `params.cnvkit_target_avg_size` added (Phase 4, second cut, 2026-08-31, see "First run — findings" below): `null` on `full`, `10000000` (10Mb bins) on `dev`.

---

## First run — findings (2026-08-31)

**Real bug, not a container/flag mistake this time: `CNVKIT_BATCH` exited status 0 but never wrote a `.cns` file**, so Nextflow correctly failed the process on `Missing output file(s) *.cns`. The container pulled and ran fine — this is not a repeat of the htslib/hap.py wrong-tag pattern.

The actual log trail: while building its internal reference from the normal BAM (`COLO829BL.dedup.bam`), CNVkit logged `Targets: 58496 (100.0%) bins failed filters (log2 < -5.0, log2 > 5.0, spread > 1.0)` — every single auto-computed bin failed CNVkit's own reference-quality QC. It still wrote `reference.cnn` anyway (with all 58,496 bins marked bad). Then, while `fix`-ing the tumour sample against that reference, it logged `Keeping 0 of 58496 bins` and `WARNING: most bins have no or very low coverage; check that the right BED file was used`, wrote a `.cnr` with 0 regions, and `segment` — given literally nothing to segment — produced no `.cns` output at all, without raising an error Nextflow could see (hence exit 0, missing file).

Root cause, confirmed via CNVkit's own `fix.py` source (not guessed): the `"Keeping %d of %d bins"` line comes from `mask_bad_bins()`, an **always-on** filter against the reference's own per-bin log2/spread statistics — it has nothing to do with `--drop-low-coverage` (that flag only ever gets a chance to act on the tumour side, during `segment`, and here there was nothing left by then for it to matter). The real problem is upstream: CNVkit's `autobin` computed ~58,496 bins (~53kb each) for the whole GRCh38 genome — a reasonable size at realistic WGS depth — but the dev-profile subsample only has ~20,000 total mapped reads spread across the *entire* 3.1Gb genome (confirmed from the log's own `Percent reads in regions: 81.104 (of 20024 mapped)` and `docs/PHASE1_NOTES.md`'s documented "10,000-read-pair subsample taken randomly across the whole genome") — roughly **0.001x effective genome-wide coverage**. At 53kb bins, the expected read count per bin is a small fraction of one, so essentially every bin's log2 ratio and spread are pure sampling noise, and 100% legitimately fail QC. This is a genuine statistical floor, not a bug in this pipeline's code or a wrong flag — the same shape of finding as Phase 1's "duplication rate needs real sequencing depth to mean anything" and Phase 2's "unsharded genome-wide Mutect2 doesn't fit this machine, at any memory size."

**Fix:** added `params.cnvkit_target_avg_size` (`nextflow.config`), set to `10000000` (10Mb) on the `dev` profile only, `null` (let CNVkit's own `autobin` decide, as it's designed to at real depth) on `full`. `modules/cnvkit.nf`'s `CNVKIT_BATCH` now takes this as an input and conditionally adds `--target-avg-size` to the command. 10Mb bins over the ~3.1Gb genome means roughly 300 bins instead of ~58,500 — concentrating the same ~20,000 reads into far fewer, far larger bins, so each one has an expected read count in the tens rather than a fraction of one, which should plausibly clear the reference-quality filter. This is a `dev`-profile-only rescue for this specific, unusually sparse smoke-test data, not a general CNVkit tuning change — it's deliberately left unset on `full`, where real WGS depth shouldn't need it. Not yet re-run to confirm the fix actually works — see "Real, unverified risks" below, this is still an open question, just a better-reasoned attempt than the original default.

---

## Design choices made, and why (so you can sanity-check them)

1. **`-m wgs` (whole-genome mode), no target/antitarget BED.** Confirmed via CNVkit's own docs (`nonhybrid.html`): in WGS mode, the reference genome's sequencing-accessible regions ("access" BED) are used as CNVkit's "targets" and are computed on the fly if not supplied — there's no separate access-BED prep step needed, unlike a capture/exome workflow.
2. **No auto-detect-or-build wrapper, unlike `.fai`/`.dict`/the truth VCF.** `INDEX_FASTA`, `CREATE_SEQUENCE_DICTIONARY`, and `PREPARE_TRUTH_VCF` all follow an "auto-detect existing file, else build it" pattern because those artifacts are expensive-ish, shared, and reused identically across every run. CNVkit's own internal reference-building (from the normal BAM) and access-BED computation are comparatively cheap, per-run, and tied to whatever tumour/normal pair is being analysed that run — there's nothing worth caching outside the process, so `CNVKIT_BATCH` just does it inline every time, matching CNVkit's own designed usage pattern rather than working against it.
3. **`-y`/`--male-reference` included, verified rather than guessed.** Getting sex-chromosome ploidy handling wrong silently biases every chrX/chrY copy-number call (CNVkit defaults to assuming a female/XX reference otherwise). Checked before writing any code, not assumed: [ATCC's own COLO 829 (CRL-1974) product page](https://www.atcc.org/products/crl-1974) lists the donor as a 45-year-old male. `-y` is passed to both `CNVKIT_BATCH` (which builds the reference internally) and `CNVKIT_CALL` — CNVkit's own docs are explicit that the same flag must be repeated at `call` time for consistent handling, not just at reference-build time.
4. **`--drop-low-coverage` included.** CNVkit's own documentation recommends this specifically "for tumor samples" — it drops bins at/near zero read depth from segmentation rather than treating them as real (very negative) log2-ratio signal. This matters more than usual here: the only tumour BAM that exists so far is Phase 2's `dev`-profile result, a 10,000-read-pair subsample spread across the whole genome, so the overwhelming majority of whole-genome bins will have exactly zero coverage.
5. **No `--annotate refFlat.txt`, deferred.** CNVkit's docs describe a gene-annotation database as needed for *readable per-target gene labels* in the output, not for `batch` to run at all. Sourcing, downloading, and verifying a GRCh38 refFlat.txt is a real task (another resource to pin a version/URL for in `docs/data_sources.md`, in the same spirit as the GATK resource bundles) — deferred rather than rushed in for a first cut whose query BAM has almost no real coverage to label meaningfully anyway. Outputs will carry generic/positional bin labels instead of gene symbols until this is revisited.
6. **No `--scatter`/`--diagram`, deferred.** Both are valid, well-documented `batch` flags for QC plots. Deliberately not added this time — Phase 3 just got burned by adding `som.py --happy-stats` (a flag that turned out to need companion flags this project wasn't ready to pick correctly) without checking what it actually required to succeed first. Same risk shape here: better to prove the base `.cnr`/`.cns`/`.call.cns` path runs cleanly first, then add plotting flags once that's confirmed, rather than risk a repeat of exactly the mistake just made and documented in `docs/PHASE3_NOTES.md`.
7. **Two processes (`CNVKIT_BATCH` then `CNVKIT_CALL`), not one.** Confirmed via CNVkit's own pipeline documentation: `batch` does NOT run the `call` step automatically — it stops at segmented-but-not-yet-integer `.cns` output. `call` is a genuinely separate step, so it gets its own process, consistent with this pipeline's one-tool-invocation-per-process style elsewhere (e.g. Mutect2 → LearnReadOrientationModel → FilterMutectCalls as three processes, not one).
8. **`--target-avg-size`, dev-profile-only, added after execution.** See "First run — findings" above for the full trail. Deliberately implemented as a profile-conditional param (`params.cnvkit_target_avg_size`, mirroring the existing `interval_list`/`dev_interval_list` pattern already used for Mutect2) rather than a fixed value baked into the module — this is a rescue for this specific, unusually sparse smoke-test subsample, not a general-purpose CNVkit tuning recommendation that should also apply once real WGS-depth data exists.

---

## Real, unverified risks worth flagging before you run this

1. **Container tag — confirmed working.** `quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0` pulled and ran successfully on the first execution (2026-08-31) — unlike htslib and hap.py in Phases 2-3, this guess turned out right. No longer an open risk.
2. **Whether 10Mb bins are actually enough to clear CNVkit's reference-quality filter is unverified.** The `--target-avg-size` fix in "First run — findings" above is a reasoned, arithmetic-backed attempt (concentrating ~20,000 reads into ~300 bins instead of ~58,500), not a confirmed working value — it hasn't been re-run yet. It's plausible 10Mb bins are still too small at this subsample's ~0.001x effective coverage and still fail the same filter, in which case the honest next move is larger bins still (or accepting, the same way Phase 2 did for unsharded Mutect2, that this specific validation needs real-depth data and isn't achievable from this subsample at any bin size) — not repeatedly guessing bigger numbers. Report back whatever the next run actually shows.
3. **Memory/CPU not preemptively tuned.** CNVkit inherits the flat per-process default (4 cpus / 8GB from the active profile). Given this project has now hit real GATK-in-Docker and bwa-mem2 memory surprises more than once (Phases 1-2), this is flagged rather than assumed fine — a `dev`-profile run against a near-empty BAM (even with far fewer, larger bins now) is expected to be lightweight, so this risk mostly matters once a real `full`-profile, real-coverage BAM exists.
4. **Runtime of the on-the-fly access-BED computation across the whole GRCh38 reference is unmeasured.** CNVkit computes this by scanning the reference FASTA for N-gaps every time `batch` runs (no caching, per design choice #2 above) — the first run's log shows this completing (all 24 chromosomes' sequences extracted) in a reasonable time alongside the rest of the process, so this turned out not to be a practical concern, though it wasn't independently timed in isolation.

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
    --truth_set_vcf         truth_set/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf \
    -resume
```

`-resume` reuses everything already computed through Modules 1-5 (Phase 3 left that fully cached and signed off) — this should only actually re-execute `CNVKIT_BATCH` (its inputs changed — new `--target-avg-size` arg — so Nextflow's cache correctly treats it as needing a rerun) and `CNVKIT_CALL`.

Expected new outputs:
- `results/cnvkit/*.cnr` — bin-level log2 copy ratios (should now have ~300 regions, not 0)
- `results/cnvkit/*.cns` — segmented (but not yet integer) copy-number regions
- `results/cnvkit/*.call.cns` — integer copy-number calls (the clinically-meaningful output)
- `results/cnvkit/*.cnn` — coverage files CNVkit uses internally to build its per-run reference

**What a successful run here would prove:** CNVkit's WGS tumour/normal batch pipeline executes against real BAMs and a real reference without erroring even at this subsample's very low effective coverage, and `call` correctly consumes `batch`'s segment output. **What it would NOT prove:** anything about COLO829's actual copy-number profile — same "code correctness, not real signal" split established in every phase so far, since the only tumour BAM that exists is the near-empty `dev`-profile subsample, now just binned much more coarsely. Expect `.call.cns` to show mostly (or entirely) neutral/diploid calls, not real amplifications or deletions.

Don't be surprised (or alarmed) by a flat, uninteresting `.call.cns` — check that the pipeline **ran to completion without errors** first, same split every phase so far has established. If `CNVKIT_BATCH` still fails on `Keeping 0 of N bins`/`Missing output file(s) *.cns` even at 10Mb bins, that's a real, more significant finding (see "Real, unverified risks" #2 above) — paste the error back and we'll diagnose it the same way every other real bug in this project has been diagnosed: from the actual failure, not a guess.
