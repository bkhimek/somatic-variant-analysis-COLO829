# Mutect2 interval scatter/gather: design notes

**Date:** 2026-09-01
**Status:** Designed and implemented against GATK's own documentation, Broad's production `mutect2.wdl`, and nf-core/sarek's actual Nextflow source (all fetched and read directly, not guessed or reasoned from memory). **Not yet run.** Same "review isn't the same as execution" caveat as every other phase/module in this project — this doc records the research trail and the choices it justifies; it does not claim the code works until it's actually been executed.

**Why this exists:** Phase 2 (`docs/PHASE2_NOTES.md`) found that a single, unsharded Mutect2 invocation against the whole genome doesn't fit this machine at any memory size tried, and deferred a real fix — interval scatter/gather — until real full-genome data actually made it necessary. Phases 3 and 4 both signed off explicitly on the basis that their results were "code correctness, not real signal," deferring meaningful precision/recall and copy-number findings to the same trigger. That point has been reached: real full-genome tumour/normal BAMs already exist from Phase 1's real full-genome alignment run. This module makes it possible to actually run Mutect2 against them.

---

## Research trail (sourced, not guessed)

**1. GATK's own recommended interval-splitting approach.**
- `SplitIntervals` is the scatter tool (`--scatter-count`/`-scatter`, `--subdivision-mode`/`-mode`, `-L`, `-O <output directory>`). Confirmed via GATK's own doc page. `--subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION` keeps each input interval intact instead of splitting inside one — GATK's own doc recommends this for assembly-based callers like Mutect2, to avoid an analytical edge artifact at a shard boundary that falls mid-interval.
- Given no `-L` at all, `SplitIntervals` splits the whole reference by raw base count, which can cross an assembly-gap (N) region arbitrarily (confirmed via a GATK developer's own forum answer). Broad's own practice — and how their published `wgs_calling_regions.hg38.interval_list` was actually built — is to first run `ScatterIntervalsByNs` (`-R`, `-O`, `-OT`/`--OUTPUT_TYPE ACGT`) to get an N-gap-aware, ACGT-only interval list, then feed that into `SplitIntervals`.

**2. Broad's production `mutect2.wdl`** (verified directly against `raw.githubusercontent.com/gatk-workflows/gatk4-somatic-snvs-indels/master/mutect2.wdl`):
- `SplitIntervals` scatter task calls `gatk SplitIntervals -R ref -L intervals -scatter ${scatter_count} -O interval-files` — `scatter_count` is a required workflow input, not a hardcoded default.
- VCF gather uses **`MergeVcfs`**, explicitly not `GatherVcfs` — the WDL's own inline comment: *"using MergeVcfs instead of GatherVcfs so we can create indices — WARNING 2015-10-28: GatherVcfs Index creation not currently supported when gathering block compressed VCFs."*
- Stats gather uses **`MergeMutectStats -stats f1 -stats f2 ... -O merged.stats`** (repeated `-stats` flag).
- Per-shard `f1r2.tar.gz` files are **not** pre-merged as a separate step — `LearnReadOrientationModel` is called once, directly on the full array of per-shard f1r2 tarballs (`-I file1 -I file2 ...`), and produces one artifact-priors tarball. GATK's own community guidance for scattered Mutect2 runs says the same thing explicitly: *"if you are scattering Mutect2 ... you must input the --f1r2-tar-gz output from each Mutect2 scatter to LearnReadOrientationModel."*
- `FilterMutectCalls` takes the single merged VCF, single merged stats file, one contamination table, and one artifact-priors tarball — everything gathered to one-per-sample before filtering.

**3. nf-core/sarek's actual implementation** (verified against `subworkflows/local/bam_variant_calling_somatic_mutect2/main.nf` and the modules it includes): sarek does **not** use `SplitIntervals` — it has its own hand-written AWK-based interval chunker targeting a runtime budget (`params.nucleotides_per_second`). Its gather side matches Broad's WDL exactly, though: `GATK4_MERGEVCFS` and `GATK4_MERGEMUTECTSTATS`, and a single `GATK4_LEARNREADORIENTATIONMODEL` call fed every shard's f1r2 tarball via repeated `--input`. Sarek's custom scatter approach wasn't adopted here — `SplitIntervals` is GATK's own documented tool for exactly this and needs no hand-written chunking logic, and both reference pipelines agree on the gather side regardless of which scatter method they use.

**4. Merge tools before `FilterMutectCalls`, confirmed against both reference pipelines:** VCFs via `MergeVcfs`, stats via `MergeMutectStats`, f1r2 via a single multi-input `LearnReadOrientationModel` call (not merged as its own step) — this project's `modules/mutect2.nf` follows exactly this shape.

**5. Single-machine/Docker-specific gotchas:** GATK's own Mutect2 FAQ says one CPU per instance is enough — parallelism is meant to come from running multiple shards concurrently, not from multithreading one. No documented GATK-level locking/concurrency issue was found for multiple Mutect2 JVMs reading the same reference/BAM concurrently (each opens independent file handles; BAM random access via `.bai` is designed for concurrent readers). The practical constraint on this machine is host resource contention — N concurrent JVMs each need their own heap — not GATK correctness. **No GATK-documented "right" scatter count for single-machine parallelism exists** — this is a deployment choice, not a sourced parameter, and is flagged as such rather than presented as researched fact (see "Design choices" #5 below).

---

## Design choices made, and why

1. **`SplitIntervals`, not sarek's custom AWK chunker.** It's GATK's own tool, built for exactly this, needs no hand-written logic to maintain, and both Broad's WDL and (for the gather side) sarek agree on how its output should be merged back together.
2. **`ScatterIntervalsByNs` runs only when `params.interval_list` is unset (the `full` profile).** On `dev`, `params.interval_list` is already the small, curated melanoma-driver-gene BED — real gene coordinates, not raw genome scanning, so there's no N-gap risk to guard against and running an extra whole-genome N-scan there would be pure overhead. `workflows/somatic.nf` branches on whether `params.interval_list` is set, not on profile name directly, keeping this consistent with the same NO_FILE-sentinel-adjacent pattern already used for `GET_PILEUP_SUMMARIES`' interval restriction.
3. **`MergeVcfs`, not `GatherVcfs`.** Directly inherited from Broad's own stated reasoning (index creation support), not independently re-derived — no reason to make a different call than the pipeline this project already treats as the reference implementation for exactly this step.
4. **`LearnReadOrientationModel` fed a collected list, not run per-shard-then-merged.** Both reference pipelines do it this way, and it matches the tool's own documented multi-`-I` capability — there is no separate "merge the f1r2 tarballs" tool to reach for, because none is needed.
5. **`mutect2_scatter_count`: 2 on `dev`, 20 on `full`, both profile-conditional params.** The `dev` value exists purely to exercise the >1-shard merge path cheaply (proving `MergeVcfs`/`MergeMutectStats`/`LearnReadOrientationModel` actually combine more than one shard correctly) before trusting it on a real whole-genome run — not for real parallelism benefit, since `dev`'s calling region is already tiny. The `full` value (~155Mb of calling region per shard, given a ~3.1Gb genome) is a **reasoned starting guess, explicitly not a sourced GATK/Broad/sarek number** (research finding #5 above states plainly that no such universal number exists) — flagged as an open risk below, in the same spirit as this project's other resource numbers (BWA_MEM2_ALIGN's memory, CNVkit's bin-size override) that were set as a reasoned starting point and corrected by real execution rather than never guessed at all.
6. **`MUTECT2`'s `cpus` dropped to 1 (from the flat 4-cpu profile default), set directly on the process rather than via a profile `withName` override.** Backed by GATK's own FAQ guidance (one CPU per instance is enough) — freeing more of the executor's shared cpu pool (`nextflow.config`, currently capped at 4 cpus / 26GB total) for multiple shards to actually run concurrently, rather than one 4-cpu Mutect2 task blocking everything else. Set on the process itself (not per-profile) because this is a property of how Mutect2 itself parallelizes, true on both `dev` and `full`, not a machine-specific tuning knob like the BWA_MEM2_ALIGN memory override is.
7. **Exact, shard-derived output filenames on `MUTECT2` (`${interval_shard.baseName}.unfiltered.vcf.gz` etc.), not a fixed name.** Every shard task writes into its own isolated Nextflow work directory, so same-named outputs don't collide there — but `publishDir` copies every task's output into one shared results directory, where a fixed name across N shards would silently overwrite itself N-1 times. This is the same category of mistake as Phase 4's `CNVKIT_CALL` bug (`docs/PHASE4_NOTES.md`'s "Second run — findings"), just encountered from the opposite direction (collision on the way out via `publishDir`, rather than an unexpectedly-multi-file glob on the way in) — applied proactively here rather than waiting to hit it again.
8. **Contamination estimation (Module 3 / `GET_PILEUP_SUMMARIES`/`CALCULATE_CONTAMINATION`) is deliberately left unscattered.** The Phase 2 finding that was actually confirmed and deferred was specifically about unsharded genome-wide *Mutect2*, not `GetPileupSummaries` — there is no confirmed finding that the contamination step has the same problem. Scattering it too, without evidence it needs it, would be scope creep beyond what this rework set out to fix; if a real `full`-profile run shows `GET_PILEUP_SUMMARIES` failing the same way, that becomes its own documented finding at that point, the same "don't fix what hasn't been proven broken" discipline this project has followed throughout.

---

## Real, unverified risks worth flagging before you run this

1. **None of this has actually been executed yet.** Every choice above is sourced from GATK's own docs and two real reference pipelines, but this project's own established pattern (five phases running now) is that execution finds real bugs review doesn't catch — Nextflow glob/list-vs-single-file mismatches (Phase 4), JVM heap sizing under Docker (Phase 2), statistical floors at low coverage (Phase 4) — none of which were visible from source-reading alone before they were hit. Treat this design as unverified until a real run against real data confirms it.
2. **`mutect2_scatter_count = 20` on `full` is an unresearched starting guess**, explicitly not a GATK/Broad/sarek-documented number (see research finding #5). It may be far too few (if per-shard runtime is still too long) or too many (if per-shard task overhead or disk I/O contention on this machine's WSL2 virtualized disk dominates) — expect to tune this after a real attempt, the same way BWA_MEM2_ALIGN's memory setting went through several real, execution-driven corrections in Phase 1.
3. **`cpus 1` on `MUTECT2` changes concurrency, not memory** — the flat 8GB-per-task profile default is unchanged, so how many shards can actually run at once is still bounded by the executor's 26GB total pool (`nextflow.config`) as much as by the cpu count. Whether 20 full-profile shards at ~8GB each (only ~3 fitting in the pool at a time) turns out to be a reasonable balance of parallelism vs. wall-clock time is unverified.
4. **`GET_PILEUP_SUMMARIES`/`CALCULATE_CONTAMINATION` (Module 3) has never actually been run against the whole, unrestricted genome** — only against the `dev` profile's gene-panel-restricted BAMs (Phase 2's sign-off basis). Whether it hits a similar "doesn't fit this machine" wall the way unsharded Mutect2 did is genuinely unknown, not assumed fine — design choice #8 above deliberately left it unscattered because there's no confirmed finding yet that it needs to be, but a real `full`-profile run could surface exactly that as a new finding.
5. **CNVkit (Module 6) and `som.py` benchmarking (Module 5) both still run unscattered against the same real full-genome BAMs.** This rework only touches Mutect2. If a real `full`-profile run gets far enough to reach `CNVKIT_BATCH` or `SOMPY_BENCHMARK` against real full-depth coverage (rather than the near-empty `dev`-profile subsample both were signed off against), that's new territory neither module has been exercised in — worth watching for in the run output, not assumed to "just work" because the `dev`-profile version did.

---

## How to run it

This hasn't been tried yet — this section is what running it *should* look like, not a confirmation that it does.

```bash
conda activate nextflow
cd ~/projects/somatic-variant-analysis-COLO829

nextflow run main.nf -profile docker,full \
    --tumour_reads_1        <real full-genome tumour FASTQ R1> \
    --tumour_reads_2        <real full-genome tumour FASTQ R2> \
    --normal_reads_1        <real full-genome normal FASTQ R1> \
    --normal_reads_2        <real full-genome normal FASTQ R2> \
    --reference_fasta       reference/Homo_sapiens_assembly38.fasta \
    --panel_of_normals      reference/1000g_pon.hg38.vcf.gz \
    --germline_resource     reference/af-only-gnomad.hg38.vcf.gz \
    --common_biallelic_sites reference/small_exac_common_3.hg38.vcf.gz \
    --truth_set_vcf         truth_set/COLO-829-NovaSeq--COLO-829BL-NovaSeq.snv.indel.final.v6.annotated.vcf \
    -resume
```

`-resume` should reuse Phase 1's already-completed full-genome `BWA_MEM2_ALIGN`/`MARK_DUPLICATES` output if it's still on disk — worth confirming those dedup BAMs still exist before this run, since Phase 1 signed off in a previous session and this project has already had one instance this session (Phase 4's third run) of a stale local file/state assumption costing a wasted round trip.

Before running the real `full` profile, it's worth first running `-profile docker,dev` once to confirm the scatter/gather DAG itself is wired correctly (2 shards, small/fast) — same "prove the DAG cheaply before the expensive real run" discipline as every phase before this one — since `dev`'s `mutect2_scatter_count = 2` exists specifically for that purpose (design choice #5 above).

**What a successful `full` run would prove:** that real interval scatter/gather lets Mutect2 actually complete against the real full-genome BAMs on this machine, closing out the risk Phase 2 deferred, and unlocking real precision/recall numbers for Phase 3's benchmarking and real depth-driven CNV calls for Phase 4's CNVkit output (both currently only proven as "wiring correct," not "biologically meaningful"). **What it would NOT automatically prove:** that `mutect2_scatter_count = 20` is the *right* number, or that Modules 3/5/6 handle real full-depth data cleanly — see "Real, unverified risks" above.
