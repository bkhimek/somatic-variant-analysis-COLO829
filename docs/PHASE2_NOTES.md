# Phase 2 — Contamination Estimation and Mutect2 Somatic Calling (Modules 3–4): implementation notes

**Date:** 2026-08-30
**Status:** Ran successfully on the `dev` profile (driver-gene-restricted) first try; the `full`-profile (unrestricted) rerun has so far surfaced two real bugs, both fixed — see "First real run — findings" and "Second `full`-profile attempt — Mutect2 JVM heap OOM" below. Originally built and reviewed in this session, following the same GATK Best Practices tumour-normal design the repo's existing scaffolding already signposted (`panel_of_normals`, `germline_resource`, `common_biallelic_sites` params, and the `qc.contamination` warn/fail thresholds in `conf/resources.config`, were all already declared in Phase 0/1, unused until now). Reviewed carefully for the DSL2 gotchas this project has already been burned by twice (Phase 1's channel-cardinality bugs, the Nextflow 26.x directive-closure issue) — but as expected, review alone didn't catch everything; each `full`-profile rerun so far has found one more real bug the `dev`-profile run's tiny restricted scope had masked.

---

## First real run — findings (2026-08-30)

**`dev`-profile run (driver-gene-BED-restricted) succeeded first try** — all 9 processes (`INDEX_VCF`, `GET_PILEUP_SUMMARIES`, `CALCULATE_CONTAMINATION`, `MUTECT2`, `LEARN_READ_ORIENTATION_MODEL`, `FILTER_MUTECT_CALLS`, plus cached Modules 1-2) completed without error. Output, as expected per the "Interval BED" caveat above: `contamination.table` showed `contamination=0.0, error=1.0` (an `error` of 1.0 signals "no reliable estimate," not "genuinely zero contamination" — there simply weren't enough informative sites in those 8 tiny gene windows given the sparse subsample), `segments.table` had zero segments, and `filtered.vcf.gz` had zero variant calls. This confirms the DAG/containers/GATK command lines are wired correctly — not that there's real signal to see yet, exactly the split Phase 1 already established with its chr21 test.

**Real bug found via the `full`-profile (unrestricted) rerun: `GetPileupSummaries` requires `-L`/`--intervals` — it is not optional.** The original design treated the interval restriction as optional for `GET_PILEUP_SUMMARIES` (same NO_FILE-sentinel pattern that correctly makes `-L` optional for `MUTECT2`), reasoning that omitting it entirely would mean "scan everywhere" on the `full` profile. That's valid for Mutect2, but GATK 4.5.0.0's `GetPileupSummaries` hard-fails without it: `A USER ERROR has occurred: Argument intervals was missing: Argument 'intervals' is required`. The `dev`-profile run never hit this because it always passed the driver-gene BED as `-L`, satisfying the requirement — the bug was only reachable through the `full`-profile path, which is exactly why trying the unrestricted rerun (rather than stopping at the dev-profile "ran clean" result) was worth doing.

**Fix:** standard GATK Best Practices pattern — always pass the common-biallelic-sites VCF itself as the baseline `-L` (those are the only positions `GetPileupSummaries` ever needs to check pileups at regardless of any further restriction), and additionally intersect with a real `interval_list` when one is set (`-L <bed> -L <common_sites_vcf> --interval-set-rule INTERSECTION`) rather than replacing the sites restriction outright. This satisfies GATK's requirement unconditionally and is the intended, efficient way to run this tool (traversing the whole genome checking every possible position against a sites VCF would be wasteful even if it were technically allowed).

---

## Second `full`-profile attempt — Mutect2 JVM heap OOM (2026-08-30)

With the `GetPileupSummaries` fix landed, the unrestricted `full`-profile rerun got further: `GET_PILEUP_SUMMARIES` completed successfully for both samples (2 of 2). `MUTECT2` then failed after 5.82 minutes with exit code 137:

```
[...] Mutect2 done. Elapsed time: 5.82 minutes.
Runtime.totalMemory()=2147483648
java.lang.OutOfMemoryError: Java heap space
      at htsjdk.samtools.SAMTextHeaderCodec$ParsedHeaderLine.<init>(...)
      ...
```

**Diagnosis:** `Runtime.totalMemory()=2147483648` is exactly 2GiB — the JVM had only self-allocated a 2GB heap despite `MUTECT2` running in this profile's flat 8GB-per-process container allocation (no per-process `memory` override existed for any Module 3–4 process). This is a well-known, well-documented property of running `gatk` under Docker: the wrapper script's own default heap-sizing logic doesn't reliably detect/scale to the actual container memory limit, so it falls back to a JVM default (often a fraction of *host* memory as seen by `free`, or a hardcoded default) that has no relationship to what the container was actually given. The crash itself is direct evidence of this — not a repeat of the earlier "~19GB to use a bwa-mem2 index" mistake, which was an unverified number pulled from an unreliable memory rather than something a crash log confirmed.

**Fix:** pass an explicit heap size via `gatk --java-options "-Xmx6g" <tool>` — this is GATK's own documented mechanism for exactly this problem. `6g` leaves ~2GB of headroom below the 8GB container allocation for off-heap/native memory, JVM overhead (metaspace, thread stacks, JIT code cache), and HTSJDK's own native buffers, which is standard practice for GATK-in-Docker sizing (heap ≈ 75% of container memory, not 100%).

Since this is a systemic property of the `gatk` wrapper under Docker — not something specific to Mutect2's workload — the same `--java-options "-Xmx6g"` flag was applied proactively to **every** GATK invocation in the pipeline, not just `MUTECT2`: `LEARN_READ_ORIENTATION_MODEL` and `FILTER_MUTECT_CALLS` (`modules/mutect2.nf`), `GET_PILEUP_SUMMARIES` and `CALCULATE_CONTAMINATION` (`modules/contamination.nf`), and `CREATE_SEQUENCE_DICTIONARY` and `INDEX_VCF` (`modules/reference_prep.nf`). `GET_PILEUP_SUMMARIES` had already completed successfully without it on this run (its workload is lighter and apparently fit inside whatever heap the JVM defaulted to), and the small resource-VCF indexing/dict-creation steps are unlikely to ever need 6GB — but leaving them unfixed would just mean the same failure loop resurfacing on a different process on some future, larger run, for no benefit (a fixed 6GB heap request costs nothing when the actual workload needs less).

**Not yet re-verified:** this fix has not yet been rerun. `CALCULATE_CONTAMINATION` was still "0 of 1" when this run aborted, so its outcome on the unrestricted `full`-profile data remains unknown until the rerun completes.

---

## What's implemented

- `modules/reference_prep.nf` — `INDEX_FASTA` (samtools faidx), `CREATE_SEQUENCE_DICTIONARY` (gatk), `INDEX_VCF` (generic tabix indexer, reused for all three resource VCFs)
- `modules/contamination.nf` — `GET_PILEUP_SUMMARIES` (per sample), `CALCULATE_CONTAMINATION` (tumour + matched normal together)
- `modules/mutect2.nf` — `MUTECT2` (tumour + normal together), `LEARN_READ_ORIENTATION_MODEL`, `FILTER_MUTECT_CALLS`
- `workflows/somatic.nf` — extended (not forked) to wire all of the above into the same `SOMATIC` workflow, per the established "Phase 2+ extends the same workflow" design
- `data/gene_lists/dev_intervals.bed` — real GRCh38 coordinates for the 8 melanoma driver genes (see "Interval BED" below) — this was an open TODO since Phase 0/1, now resolved
- `assets/NO_FILE` — sentinel placeholder enabling an "optional" interval-list input (nf-core's standard pattern) so the same processes work whether `--interval_list` is set (dev) or null (full)
- `main.nf` — `panel_of_normals`, `germline_resource`, `common_biallelic_sites` added to the required-params fail-fast check; all three now passed through to `SOMATIC`
- `nextflow.config` — **bug found and fixed while wiring this up:** `conf/resources.config` (created in Phase 1, holding the contamination warn/fail thresholds) was never actually `includeConfig`'d anywhere, so `params.qc.contamination.*` would have been undefined the first time `CALCULATE_CONTAMINATION` ran. Added the missing `includeConfig 'conf/resources.config'` line.

---

## Reference/resource prep — auto-detect-or-build, extending the Phase 1 pattern

GATK's Mutect2/GetPileupSummaries/FilterMutectCalls all need a `.fai` and `.dict` alongside the reference FASTA, and a `.tbi` alongside every resource VCF — flagged as an open TODO since `docs/PHASE0_FINDINGS.md` §7 action 8 and never actually built. `workflows/somatic.nf` now auto-detects and builds whichever of these is missing, mirroring the bwa-mem2 pre-built-index pattern from Phase 1 — **but these are all cheap (seconds), not a repeat of that memory saga.**

One deliberate simplification: unlike the `.fai`/`.dict`/bwa-mem2-index checks (which each skip an expensive rebuild if a file already exists), the three resource VCFs are always (re-)indexed via one `INDEX_VCF` call over a combined channel, rather than auto-detecting each one individually. Reason: Nextflow only cleanly supports invoking the same process multiple times in one workflow scope if each call's result is captured separately rather than read back through the shared, ambiguous `PROCESS.out` after more than one invocation — the auto-detect-per-VCF version would have needed three separate `INDEX_VCF` calls, hitting exactly that ambiguity. Since `gatk IndexFeatureFile` takes seconds (not bwa-mem2's problem), it wasn't worth the extra complexity to preserve a "skip if already indexed" optimization here — `-resume` already skips it on every run after the first regardless.

---

## Interval BED — real coordinates, but read the caveat below before trusting the numbers

`data/gene_lists/dev_intervals.bed` now has real GRCh38 coordinates for all 8 genes in `melanoma_genes.tsv` (BRAF, NRAS, CDKN2A, PTEN, TP53, TERT, NF1, KIT), resolved via Ensembl's REST API (`https://rest.ensembl.org/lookup/symbol/homo_sapiens/<GENE>`) rather than guessed — each lookup's `assembly_name` was confirmed `"GRCh38"` before use, and coordinates were padded ±2000bp per `melanoma_genes.tsv`'s own TODO note ("+/- flanking region"). Chromosome names are `chr`-prefixed to match `Homo_sapiens_assembly38.fasta`'s convention (Ensembl's own API returns bare numbers).

**TERT note:** TERT is minus-strand, so its well-known recurrent promoter mutations (C228T/C250T) sit just upstream of the gene's *high*-coordinate end, not the low-coordinate start. Rather than assert an exact hotspot coordinate from memory (unverified), the uniform 2000bp padding is relied on to cover the promoter region without pinpointing it — if precise hotspot-level resolution is ever needed later (e.g. for a hotspot-specific report), verify the exact positions independently first.

**The important caveat — unlike Mutect2/GetPileupSummaries respecting `-L` at all (which alignment doesn't, see `nextflow.config`'s dev-profile comment), this restriction interacting with the *current* dev-profile test data is limited:** the dedup BAMs that exist right now (`results/alignment/*/*.dedup.bam`) were built from a 10,000-read-pair subsample taken **randomly across the whole genome** (`docs/PHASE1_NOTES.md`), not from reads specifically overlapping these 8 gene windows. There's no reason to expect meaningful (or any) coverage at these specific loci in that subsample. Restricting Mutect2/contamination to this BED against the *current* BAMs will most likely validate that the DAG, containers, and GATK command lines are wired correctly — not that real variant-calling behaviour looks right, since there's essentially nothing to call there. This mirrors Phase 1's chr21-test distinction (code correctness vs. real signal) — expect the same kind of split result, and don't read too much into an empty or near-empty `filtered.vcf.gz` from this specific test.

---

## Sample branching — genuinely different from Modules 1–2

Modules 1–2 process tumour and normal **identically and independently** (same process, run twice). Modules 3–4 don't: `CalculateContamination` needs the tumour's pileups analysed against the matched normal's, and `Mutect2` needs both BAMs in one invocation with `-normal <normal_sample_id>`. `workflows/somatic.nf` uses `.branch{ tumour: it[0] == 'COLO829'; normal: it[0] == 'COLO829BL' }` twice (once on the dedup BAMs, once on the pileup tables) to split the shared per-sample channel back into named tumour/normal channels for the steps that need pairing. This is exactly why `main.nf`'s hardcoded sample IDs (`COLO829`/`COLO829BL`) mattered from the start — `PHASE1_NOTES.md` flagged keeping them consistent "since Mutect2's tumour/normal pairing... key off them," and this is that pairing actually landing.

---

## Required params — Phase 1-only runs now need three more

Since Phase 2 extends the same `SOMATIC` workflow rather than making Modules 3–4 optional, `main.nf`'s required-params check now also requires `--panel_of_normals`, `--germline_resource`, and `--common_biallelic_sites` (provenance/download commands: `docs/data_sources.md` §4). A command that worked for Phase 1 alone will now fail fast with a clear missing-param message rather than silently skipping the new modules — that's intentional, not a regression.

---

## Things you need to check/fix before this actually runs

1. **This has never been executed.** Expect at least one round of real bugs, same as every phase so far — DSL2 channel-cardinality mistakes are exactly the kind of thing that only shows up at runtime (Phase 1 had several).
2. **GATK JVM heap is now explicitly sized (`-Xmx6g` on every GATK invocation), container memory is not.** All Module 3–4 processes still inherit the profile's flat 8GB container allocation — see "Second `full`-profile attempt" above for why the *heap* needed an explicit size regardless of container memory (a Docker-specific gatk-wrapper quirk, confirmed by a real crash, not a guess). If a future run OOMs at the container level (not the JVM-heap level) on real full-depth WGS data, that's a separate, legitimate reason to raise the container `memory` directive too — raise `-Xmx` alongside it if so, keeping ~25% headroom.
3. **Download the three GATK resource VCFs** (`docs/data_sources.md` §4) if you haven't already — `panel_of_normals` (`1000g_pon.hg38.vcf.gz`), `germline_resource` (`af-only-gnomad.hg38.vcf.gz`), `common_biallelic_sites` (`small_exac_common_3.hg38.vcf.gz`). `.tbi` companion index files are commonly published alongside GATK bundle VCFs at the same path + `.tbi`, but this hasn't been independently confirmed for this specific bucket — it doesn't matter either way, since `INDEX_VCF` builds them if missing regardless (see above).
4. **The contamination warn/fail thresholds (2%/5%) are explicitly provisional** — `conf/resources.config` says so directly, pending inspection of the real COLO829 contamination table. `CALCULATE_CONTAMINATION` logs a WARN/FAIL/OK line against these thresholds but doesn't fail the pipeline on it — revisit the actual numbers once you have real output.
5. **`gnomAD` version embedded in `af-only-gnomad.hg38.vcf.gz` is still unconfirmed** (`docs/data_sources.md` §4/§9, an open TODO since Phase 0) — check the VCF header once downloaded and record it, needed for `run_manifest.json`'s `gnomad_version` field eventually.

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
    -resume
```

(Paths above assume you download the three resource VCFs into the same `reference/` directory as the FASTA — adjust if you put them elsewhere.) `-resume` will reuse everything already computed for Modules 1–2 from the Phase 1 run, so this should only actually execute the new Module 3–4 steps plus the new `.fai`/`.dict`/VCF-index prep steps.

Expected new outputs:
- `results/contamination/<sample_id>/<sample_id>.pileups.table`
- `results/contamination/contamination.table`, `results/contamination/segments.table`
- `results/mutect2/unfiltered.vcf.gz` (+ `.tbi`, `.stats`, `f1r2.tar.gz`)
- `results/mutect2/read-orientation-model.tar.gz`
- `results/mutect2/filtered.vcf.gz` (+ `.tbi`, `.filteringStats.tsv`) — **this is Phase 3's benchmarking-against-truth-set target**

Per the caveat above, don't be surprised (or alarmed) if `filtered.vcf.gz` from this specific dev-profile run against the current subsample BAMs comes back empty or near-empty — check that the pipeline **ran to completion without errors** first, the same "code correctness before real signal" split Phase 1 established.
