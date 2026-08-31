# Phase 4 — Copy-number calling (Module 6): implementation notes

**Date:** 2026-08-31
**Status:** PHASE 4 SIGNED OFF 2026-08-31 — three runs, three real findings, all root-caused via execution rather than guessed. First run: a genuine statistical floor (`CNVKIT_BATCH` exited 0 without ever writing a `.cns`, because the dev-profile subsample is too sparse for CNVkit's default whole-genome bin size), fixed with a `dev`-profile-only bin-size override. Second run: `CNVKIT_CALL` failed on a 3-file glob collision (`batch` writes more `.cns`-suffixed files than this module first assumed), fixed with exact output filenames and separated `publishDir`s. Third run: the DAG completed cleanly end-to-end, but only after discovering that Nextflow's `-resume` doesn't re-derive a cached task's outputs against a changed output-declaration — it replays the exact file list recorded the first time that task ran, so the fix from run two silently didn't take effect until `CNVKIT_BATCH`'s own work directory was deleted to force a real re-execution. Once genuinely re-run, the whole DAG including `CNVKIT_CALL` completed successfully. The resulting `.call.cns` is exactly the "mechanically correct, no real signal" result predicted for this subsample — see "Third run — result" and "Phase 4 sign-off" below. Same "review isn't the same as execution" caveat as every phase so far applies — and this phase adds its own coda to it: execution isn't even proof of the CURRENT code once `-resume` is in the mix.

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

**Fix:** added `params.cnvkit_target_avg_size` (`nextflow.config`), set to `10000000` (10Mb) on the `dev` profile only, `null` (let CNVkit's own `autobin` decide, as it's designed to at real depth) on `full`. `modules/cnvkit.nf`'s `CNVKIT_BATCH` now takes this as an input and conditionally adds `--target-avg-size` to the command. 10Mb bins over the ~3.1Gb genome means roughly 300 bins instead of ~58,500 — concentrating the same ~20,000 reads into far fewer, far larger bins, so each one has an expected read count in the tens rather than a fraction of one, which should plausibly clear the reference-quality filter. This is a `dev`-profile-only rescue for this specific, unusually sparse smoke-test data, not a general CNVkit tuning change — it's deliberately left unset on `full`, where real WGS depth shouldn't need it.

---

## Second run — findings (2026-08-31)

**With the bin-size fix landed, `CNVKIT_BATCH` succeeded** (confirming the 10Mb override worked — a real, if indirect, confirmation that the "First run" theory was correct). **`CNVKIT_CALL` then failed**: `Command exit status: 2`, `cnvkit.py: error: unrecognized arguments: COLO829.dedup.call.cns COLO829.dedup.cns COLO829.dedup.call, COLO829.dedup].call.cns`.

What actually happened: `CNVKIT_BATCH`'s work directory contained **three** files matching this module's `path("*.cns")` glob — `COLO829.dedup.cns` (the plain segmented file this module expected), `COLO829.dedup.call.cns`, and `COLO829.dedup.bintest.cns`. Nextflow correctly staged all three as a 3-element list into `CNVKIT_CALL`, and interpolating a list into the script's `${cns_file}` placeholder writes all three filenames space-separated on the command line — `cnvkit.py call` took the first as its one expected positional argument and rejected the rest as unrecognized.

Where the extra two files came from is **not fully root-caused**, and that's stated plainly rather than papered over: checked CNVkit's current (GitHub master-branch) `cnvlib/commands.py` and `cnvlib/segmentation/__init__.py` source directly, and neither shows `batch` invoking `call` or `bintest` automatically — `batch`'s own documentation (`pipeline.html`, checked when this module was first designed) says the same. Two honest possibilities, neither confirmed: (1) the pinned `0.9.10` release's actual `batch.py` behaves differently from current master (bioconda's cnvkit recipe has already moved on to `0.9.13`, so real drift between the pinned version and the source checked is plausible), or (2) some flag combination here (`-y`, `--drop-low-coverage`, `--target-avg-size`) triggers additional output in ways not obvious from a source skim. Either way, this is a confirmed-by-execution fact, not a guess — the fix doesn't depend on knowing which explanation is right.

**Fix:** made `CNVKIT_BATCH`'s `cnr`/`cns` outputs exact filenames (`"${tumour_bam.baseName}.cnr"` / `"${tumour_bam.baseName}.cns"`, using the same input-BAM-basename convention CNVkit itself uses and this pipeline already relies on elsewhere) instead of glob patterns — `CNVKIT_CALL` now only ever receives the one plain segmented file, regardless of whatever else `batch` produces alongside it. Also split `CNVKIT_BATCH`/`CNVKIT_CALL` into separate `publishDir` subfolders (`${params.outdir}/cnvkit/batch`, `${params.outdir}/cnvkit/call`) — otherwise `CNVKIT_BATCH`'s own mystery `COLO829.dedup.call.cns` and `CNVKIT_CALL`'s genuine output of the identical name would have silently collided in the same output directory once both got published. The mystery `.call.cns`/`.bintest.cns` files are still captured (via the `all_outputs` catch-all) and published under `cnvkit/batch/` for reference, just no longer fed into `CNVKIT_CALL`.

---

## Third run — findings and result (2026-08-31, SUCCESS)

**The fix above didn't take effect on the first attempt to re-run it, and the reason is a genuinely new finding about how `-resume` works, not a mistake in the fix itself.** The first re-run attempt showed `CNVKIT_BATCH` as `cached: 1 ✔` and `CNVKIT_CALL` still received the exact same 3-file list as the second run, with the exact same "unrecognized arguments" error. This was surprising, since the module file on disk genuinely had the exact-filename fix in it by then.

Root cause, confirmed rather than guessed: `CNVKIT_BATCH`'s script and inputs (BAMs, reference, `target_avg_size`) hadn't changed between the second and third run attempts — only its `output:` block had. Nextflow's task cache key (used to decide whether `-resume` can reuse a task) is built from the process's script and inputs, not its output declarations, and a cache hit **replays the exact output file list recorded when that task last actually ran** rather than re-deriving it against whatever the output block currently says. So even with the corrected module in hand, `-resume` kept handing `CNVKIT_CALL` the stale 3-file list captured back when `CNVKIT_BATCH` still used the `path("*.cns")` glob. The fix only actually took effect once `CNVKIT_BATCH`'s specific work directory (`work/ee/f2e427...`) was deleted, forcing a genuine re-execution — visible in the log as `CNVKIT_BATCH` completing without the `cached` label this time.

This is worth calling out on its own, separate from the two module bugs above: **a change to only a process's `output:` block is invisible to `-resume` for any task whose script/inputs are unchanged**, until that task's own work directory is removed (or `-resume` is dropped entirely for a full re-run). Filed here as a process lesson for the rest of this project, not just a Phase-4-specific fix, since every phase so far has relied on `-resume` to keep iteration fast.

**Result, once genuinely re-run:** the full DAG completed end-to-end — `CNVKIT_BATCH (COLO829)` and `CNVKIT_CALL (COLO829.dedup)` (note the tag: a single filename, not a list, confirming the fix). `results/cnvkit/batch/` contains the expected `.cnr`/`.cns`/`.cnn`/`reference.cnn` plus the still-unexplained `.call.cns`/`.bintest.cns` from `batch` itself (harmless — no longer consumed downstream). `results/cnvkit/call/COLO829.dedup.call.cns` contains exactly one process's real output this time:

```
chromosome  start      end        gene  log2       ci_hi      ci_lo      cn  depth       probes  weight
chr1        2746290    12954384   -     -0.898137  -0.898137  -0.898137  1   3.04e-06    1       0.860655
chr8        135184567  145078636  -     0.0821563  0.0821563  0.0821563  2   4.65e-06    1       0.211697
chr9        128318407  138334717  -     0          0          0          2   3.98e-05    1       1
chr14       97186828   106883718  -     0.734048   0.734048   0.734048   4   2.72e-05    1       0.775327
chr19       12254490   24448980   -     -0.0456617 -0.0456617 -0.0456617 2   8.86e-06    1       0.640482
```

Read plainly, not over-interpreted: exactly 5 segments survive across the whole genome (the 10Mb-bin override plus `--drop-low-coverage` means every other bin had zero coverage and was dropped before segmentation even ran). Every surviving segment has `probes: 1` and `ci_hi == ci_lo == log2` — a single-bin "segment" has no internal variance to report a real confidence interval from, so this isn't CNVkit expressing confidence in these calls, it's CNVkit correctly reporting that it only had one data point each. `depth` values around 1e-5–1e-6 confirm the ~0.001x effective coverage figure from the "First run" finding is the right order of magnitude. The CN=1 (chr1) and CN=4 (chr14) calls are not evidence of a real deletion or amplification — they're what a single, near-zero-depth bin's noise looks like once forced through `call`'s log2-to-integer thresholds. This is precisely the "mechanically correct, no real biological signal" result predicted in this module's own header comment and in "How to run it" below, now confirmed rather than assumed.

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
9. **Exact output filenames instead of `*.cnr`/`*.cns` globs, and separate `publishDir` subfolders per process, added after execution.** See "Second run — findings" above. `batch` turned out to write more `.cns`-suffixed files than expected — an exact filename sidesteps that regardless of the exact cause, and separate publish subfolders prevent a same-named-file collision between `CNVKIT_BATCH`'s own extra output and `CNVKIT_CALL`'s genuine one.
10. **No code change for the "Third run" finding — it's a re-run/cache procedure issue, not a module bug.** See "Third run — findings and result" above: the exact-filename fix from item 9 was correct as written, it just needed `CNVKIT_BATCH`'s own work directory deleted once to force `-resume` to actually apply it. Nothing in `modules/cnvkit.nf` changed for this.

---

## Real, unverified risks worth flagging before you run this again

1. **Container tag — confirmed working.** `quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0` pulled and ran successfully across all three executions (2026-08-31) — unlike htslib and hap.py in Phases 2-3, this guess turned out right. No longer an open risk.
2. **The 10Mb bin-size fix — confirmed working end-to-end.** The third run's `.call.cns` (see "Third run — result" above) shows 5 real segments with sane, if uninteresting, values — not some other degenerate failure mode the coarser binning might have introduced. No longer an open risk.
3. **`CNVKIT_CALL`'s exact-filename fix — confirmed working.** The third run's tag (`CNVKIT_CALL (COLO829.dedup)`, a single filename rather than a list) and successful completion confirm it. No longer an open risk. The `-resume`-doesn't-see-output-block-changes gotcha that delayed confirming this is now documented in "Third run — findings" as its own lesson.
4. **Memory/CPU not preemptively tuned.** CNVkit inherits the flat per-process default (4 cpus / 8GB from the active profile). Given this project has now hit real GATK-in-Docker and bwa-mem2 memory surprises more than once (Phases 1-2), this is flagged rather than assumed fine — this `dev`-profile run against a near-empty, coarsely-binned BAM turned out lightweight in practice, so this risk mostly matters once a real `full`-profile, real-coverage BAM exists.
5. **Runtime of the on-the-fly access-BED computation across the whole GRCh38 reference is unmeasured in isolation.** CNVkit computes this by scanning the reference FASTA for N-gaps every time `batch` runs (no caching, per design choice #2 above) — all three runs' logs show this completing (all 24 chromosomes' sequences extracted) in a reasonable time alongside the rest of the process, so this turned out not to be a practical concern, though it wasn't independently timed on its own.

---

## Phase 4 sign-off (2026-08-31)

**Proven, by actually running the code (not just reviewing it):**
- CNVkit's `-m wgs` tumour/normal `batch` → `call` pipeline is wired correctly against real BAMs and a real GRCh38 reference: the container pulls and runs, `batch` builds its internal reference and computes bins without a supplied access/target BED, `--drop-low-coverage`/`-y`/`--target-avg-size` are all accepted and take effect, and `call` correctly consumes `batch`'s segmented output to produce integer copy-number calls.
- Two real, non-obvious execution bugs were found and fixed the same way every other phase's bugs were: root-caused from the actual failure (CNVkit's own `fix.py` source for the bin-size floor; the actual 3-file work-directory listing for the glob collision), not guessed.
- A third, more general lesson was surfaced and documented: `-resume` does not re-derive a cached task's outputs from a changed `output:` block, only from a changed script/inputs — worth remembering for any future module edit that touches only output declarations.
- `-y`/`--male-reference` is backed by an independently verified fact (ATCC's own CRL-1974 product page), not an assumption.

**Explicitly NOT proven, and deferred on the same basis Phases 2-3 used:**
- Anything about COLO829's real copy-number profile. The only tumour/normal BAM pair that exists is Phase 2's `dev`-profile ~10,000-read-pair subsample (~0.001x effective genome-wide coverage) — the 5-segment, mostly-noise `.call.cns` documented above is the correct and expected output of that input, not a real biological finding. Real CN calls need real full-genome-depth BAMs, i.e. the same "full profile + real compute" step Phase 2 already deferred for Mutect2.
- `--annotate`/`--scatter`/`--diagram` remain deferred (design choices #5, #6) — the base `.cnr`/`.cns`/`.call.cns` path is proven first, plotting/annotation flags are a follow-up once there's real signal worth visualising.
- Memory/CPU sizing for a real `full`-profile run (risk #4 above) is still just the flat per-process default, unstressed by real WGS-scale coverage.

Signed off on that basis — same split as every phase so far: mechanical correctness of the pipeline is proven; biological meaningfulness of the numbers is not, and isn't expected to be until real full-genome data exists.

---

## How to run it

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

`-resume` will reuse everything through Modules 1-5 and `CNVKIT_BATCH`/`CNVKIT_CALL` themselves, since this is now a confirmed, signed-off good state — nothing here should re-execute unless an input file actually changes. (If you do change `modules/cnvkit.nf` again in future, remember the "Third run" finding above: a change to only the `output:` block won't be picked up by `-resume` until that specific task's `work/<hash>` directory is deleted.)

Actual outputs, confirmed 2026-08-31 (see "Third run — result" above for the full readout):
- `results/cnvkit/batch/*.cnr` — bin-level log2 copy ratios (5 non-zero regions survive `--drop-low-coverage` at this coverage level)
- `results/cnvkit/batch/*.cns` — segmented (but not yet integer) copy-number regions
- `results/cnvkit/batch/*.cnn` — coverage files CNVkit uses internally to build its per-run reference (plus the still-unexplained `.call.cns`/`.bintest.cns` `batch` writes alongside the plain `.cns` — harmless, not consumed downstream)
- `results/cnvkit/call/*.call.cns` — integer copy-number calls: 5 segments, CN values of 1/2/2/4/2, each from a single-probe segment (near-zero confidence, near-zero depth) — the expected "mechanically correct, no real signal" result for this subsample, not a real CNV finding

**Confirmed proven:** CNVkit's WGS tumour/normal batch pipeline executes against real BAMs and a real reference without erroring even at this subsample's very low effective coverage, and `call` correctly consumes `batch`'s segment output. **Confirmed NOT meaningful:** anything about COLO829's actual copy-number profile — see "Phase 4 sign-off" above for the full basis.
