# Phase 1 — QC and Alignment (Modules 1–2): implementation notes

**Date:** 2026-08-29 (updated same day after the first real run; updated again 2026-08-30 after the full-genome index build)
**Status:** Code has now actually run — see "First real run — findings" below. Originally built in a Cowork sandbox with no Nextflow, no Docker, and no access to the real FASTQs or reference; I tried to install Nextflow there to lint/run it (`get.nextflow.io` and GitHub releases were both blocked by that sandbox's network policy) and couldn't, so it was only checked by manual read-through before delivery. Manual review caught two real DSL2 channel-cardinality bugs before delivery (both classic Nextflow gotchas — see `workflows/somatic.nf` comments at the `.flatMap` and `BWA_MEM2_INDEX.out.index.first()` lines), but as flagged at the time, review alone isn't the same guarantee as an actual run — and indeed, the first real run surfaced two things review missed (below).

---

## Full-genome bwa-mem2 index — built 2026-08-30, pipeline updated to use it

Built on a one-off AWS `r6i.4xlarge` instance (128GB RAM, ~$1-2 total cost, eu-north-1) rather than locally — the 32GB local machine can't build it (§ above / `docs/data_sources.md` §2: needs ~87-110GB RAM), but only needs ~19GB to *use* an already-built index, which is well within reach locally. Index build itself took **22.4 minutes** on the AWS instance (faster than the ~50min benchmark quoted earlier, likely due to more vCPUs), producing five files (`Homo_sapiens_assembly38.fasta.{0123,amb,ann,bwt.2bit.64,pac}`, ~16.9GB total) alongside the reference FASTA in `reference/`.

**`workflows/somatic.nf` updated to auto-detect this.** Rather than adding a new param to remember to set, the workflow now checks whether `${reference_fasta}.bwt.2bit.64` already exists next to whatever `--reference_fasta` points at: if it does (as is now true for the full-genome reference), `BWA_MEM2_INDEX` is skipped entirely and the existing index files are fed straight to `BWA_MEM2_ALIGN`; if not (e.g. the chr21 dev-profile slice, which has no pre-built index), it falls back to building the index in-pipeline as before. This means the same `nextflow run` command works whether or not a pre-built index exists — no extra flag needed.

**Next step, not yet done:** rerun the same 10,000-read-pair subsample (`fastq_dev_sample/`) against the *full* genome now that this index exists, to get the real mapped-fraction/duplication-rate numbers the chr21-only run couldn't provide (see "What this run does NOT validate" below) — that's what actually closes out the plan's Phase 1 exit criterion.

---

## First real run — findings (2026-08-29)

Run on a chr21-only reference slice (full-genome `bwa-mem2 index` needs ~87-110GB RAM to build — confirmed via bwa-mem2's own documentation and a matching GitHub issue report of ~80GB/50min for this exact reference; see `docs/data_sources.md` §2 note) with a 10,000-read-pair subsample of each real FASTQ. Two things surfaced that manual review didn't catch:

1. **`samtools:1.21--h50ea8bc_0` container tag — confirmed real**, resolves item 1 below; `SAMTOOLS_SORT` ran successfully.
2. **MultiQC bug, found and fixed:** `MULTIQC`'s output declared `path("multiqc_data")`, but MultiQC actually names its data directory from the `--filename` argument's stem — `--filename multiqc_report.html` produces `multiqc_report_data`, not a fixed `multiqc_data`. Nextflow failed with "Missing output file(s) `multiqc_data`" even though MultiQC itself exited 0 and worked correctly — the module just declared the wrong output path. Fixed in `modules/fastqc.nf` to `path("multiqc_report_data")`.
3. **Benign warning, no fix needed:** `WARN: The operator 'first' is useless when applied to a value channel...` on the `BWA_MEM2_INDEX.out.index.first()` line. In this pipeline's actual wiring, `reference_fasta` is built with a plain `file(params.reference_fasta)` in `main.nf` (not a multi-item queue channel), so it — and everything downstream of a process that only consumes it — is already a value channel by construction; `.first()` doesn't change behavior here. It was added as defensive best practice against the classic gotcha where a process fed by an *actual* queue channel with more than one item would otherwise "run dry" after the first item; that risk doesn't materialize in this specific wiring, but the line is harmless and left in place rather than removed for one less non-blocking warning.

Both samples aligned successfully (`BWA_MEM2_ALIGN [100%] 2 of 2`), confirming the FASTQC → BWA-MEM2 chain works end to end. Not yet confirmed at that point: `SAMTOOLS_SORT` → `MARK_DUPLICATES` output correctness (run was still in progress at the MultiQC failure) — rerun with `-resume` after the fix above to pick up from there without redoing FASTQC/alignment.

**Rerun completed successfully** (`-resume`, 2026-08-29, duration 1h28m38s — almost entirely the one-time `BWA_MEM2_INDEX` build, everything else was cache-hit or seconds): every process in Modules 1–2 now runs to completion and produces the expected files (`*.sorted.bam`, `*.dedup.bam` + `.bai`, `*.dedup.metrics.txt`, `multiqc_report.html`). **This is the code validated end to end for the first time** — the DAG wiring, container tags, and file-naming conventions are all confirmed correct.

**What this run does NOT validate — the actual `PERCENT_DUPLICATION`/coverage numbers.** This was deliberately run against a chr21-only reference slice (§ above — full-genome indexing needs ~87-110GB RAM this machine doesn't have) using reads subsampled from the real, whole-genome FASTQs. The observed metrics:

| Sample | Mapped fraction | `PERCENT_DUPLICATION` |
|---|---|---|
| COLO829 (tumour) | ~8,636 / 20,284 records (~42.6%) | 44.4% |
| COLO829BL (normal) | ~6,536 / 20,320 records (~32.2%) | 30.8% |

Both numbers look alarming in isolation (a real WGS library at 40-90% duplication, or 30-40% unmapped, would be a red flag) but **are artifacts of the test setup, not real signal**, for two compounding reasons:
1. **Mapped fraction is inflated:** these reads' true genomic origin is anywhere across the ~3.1 Gbp genome, but they can only be reported as mapped or unmapped against chr21's ~47 Mbp (~1.5% of the genome). A meaningfully higher-than-1.5%-by-chance mapped fraction is expected because BWA-MEM2 does local/soft-clipped alignment and will place a read against a partial or repeat-driven match (e.g. Alu elements, present in huge copy number genome-wide including on chr21) rather than report it unmapped outright.
2. **Duplication rate is inflated for a purely combinatorial reason:** among the reads that *do* map, they're being packed into a reference ~66x smaller (47 Mbp vs. 3.1 Gbp) than the space they'd actually occupy in a real run. Two unrelated read pairs are consequently far more likely to coincidentally start at the same position on chr21 than they would be across the whole genome — Picard's `MarkDuplicates` correctly flags that as a duplicate by its positional definition, but it doesn't reflect real PCR/optical over-amplification.

**Net: the plan's Phase 1 exit criterion — "verify duplicate rates and coverage depth are appropriate for WGS somatic calling" — is still open**, not satisfied by this run. The cheap way to close it: once a full-genome bwa-mem2 index exists (planned via a one-time AWS high-RAM instance build, `docs/data_sources.md` §2), rerun this exact same subsampled FASTQ pair against the *full* genome instead of chr21 — same fast smoke-test cost, but the resulting mapped-fraction and duplication-rate numbers would actually mean something and could be recorded here as the real Phase 1 sign-off.

## What's implemented

- `modules/fastqc.nf` — `FASTQC` (per-sample) + `MULTIQC` (aggregates tumour + normal into one report — this cross-sample aggregation is Phase 1's one explicit addition beyond reusing Project 4's QC step, per plan §11)
- `modules/alignment.nf` — `BWA_MEM2_INDEX`, `BWA_MEM2_ALIGN`, `SAMTOOLS_SORT`, `MARK_DUPLICATES`, chained per plan §6 Module 2 ("BWA-MEM2 → samtools sort → GATK MarkDuplicates")
- `workflows/somatic.nf` — wires both modules together; runs the whole chain once per sample (tumour, normal), matching the plan's explicit note that Module 2 "runs twice"
- `main.nf` — entry point; builds the two-sample channel from CLI params, fails fast with a clear message if required params are missing

Not implemented here, deliberately deferred to their own phases per the plan's build order: contamination estimation (Module 3, Phase 2), Mutect2 (Module 4, Phase 2), everything after.

---

## Things you need to check/fix before this actually runs

**1. ~~`samtools` container tag is unconfirmed~~ — CONFIRMED 2026-08-29** via both a standalone `docker pull` and the actual first pipeline run (`SAMTOOLS_SORT` executed successfully). `modules/alignment.nf`'s `SAMTOOLS_SORT` process uses `quay.io/biocontainers/samtools:1.21--h50ea8bc_0`. (Original concern: I could not browse quay.io's tag list from the build sandbox, blocked by robots.txt, same limitation hit for CNVkit/hap.py in Phase 0 — no longer relevant now it's confirmed working.)
```bash
docker pull quay.io/biocontainers/samtools:1.21--h50ea8bc_0
```
If that fails, check https://quay.io/repository/samtools/samtools?tab=tags in your browser for the real current tag and update the `container` line in `modules/alignment.nf`.

**2. No combined bwa-mem2+samtools container is used, by design.** I looked at how nf-core/modules does this (they use a Seqera Wave-hosted mulled container with a hash-only name) and deliberately didn't copy that pattern — it doesn't match this pipeline's existing convention of one plain `quay.io/biocontainers/<tool>:<version>` per tool everywhere else. Instead, `BWA_MEM2_ALIGN` writes a `.sam` file to disk and a separate `SAMTOOLS_SORT` process (different container) sorts it. Costs one extra intermediate file per sample; keeps every container single-purpose. Reasonable to revisit later if disk I/O for the intermediate SAM becomes a real bottleneck at full-WGS scale.

**3. The `dev` profile does not make alignment itself fast.** This is worth understanding before you run `-profile dev` expecting a quick smoke test. `params.interval_list` only restricts steps that take a GATK `-L` argument (Mutect2 and friends, arriving in Phase 2) — it has no effect on `BWA_MEM2_ALIGN`, which has no genomic coordinates to restrict on before alignment exists. If you want a fast Phase 1 smoke test, point `tumour_reads_1/2` and `normal_reads_1/2` at a small subsampled FASTQ pair yourself (e.g. `seqtk sample` a few thousand reads, or extract reads from a chr21 region of an existing BAM if you have one) — the pipeline won't do this for you. Full explanation is in the comment block above the `dev` profile in `nextflow.config`.

**4. BWA-MEM2 indexing runs inside the pipeline, not as a separate manual step.** `BWA_MEM2_INDEX` builds the index from whatever `params.reference_fasta` points at, every run (Nextflow's `-resume` cache will skip re-running it once it's succeeded once for that exact reference file). This means the first run against the full GRCh38 reference will spend real time and disk on indexing before alignment even starts — budget for that in addition to the alignment/sort/dedup time itself.

**5. Reference `.fai`/`.dict` are not generated yet.** Only the BWA-MEM2 index is built in Phase 1, since that's all Module 2 needs. GATK steps in Phase 2 (Mutect2, GetPileupSummaries) will need a `.fai` and a sequence dictionary (`.dict`) — add a `samtools faidx` + `gatk CreateSequenceDictionary` step when Phase 2 starts, don't assume it's already there.

**6. Sample naming is hardcoded in `main.nf`**: `COLO829` (tumour) / `COLO829BL` (normal). This matches the naming already used throughout `docs/` and `docs/run_manifest.example.json` — keep it consistent going into Phase 2, since Mutect2's tumour/normal pairing and the run manifest's `samples` block both key off these exact strings.

---

## How to run it (once you've checked the items above)

Minimal smoke test — point every required param at real files on your machine:

```bash
conda activate nextflow
cd ~/projects/somatic-variant-analysis-COLO829

nextflow run main.nf -profile docker,dev \
    --tumour_reads_1  /path/to/COLO829_R1.fastq.gz \
    --tumour_reads_2  /path/to/COLO829_R2.fastq.gz \
    --normal_reads_1  /path/to/COLO829BL_R1.fastq.gz \
    --normal_reads_2  /path/to/COLO829BL_R2.fastq.gz \
    --reference_fasta /path/to/GRCh38.fasta \
    -resume
```

If it fails on the `SAMTOOLS_SORT` container pull, that's item 1 above — fix the tag and rerun with `-resume` (BWA-MEM2 indexing and alignment results already computed will be reused, not redone).

Expected outputs once it completes:
- `results/qc/multiqc/multiqc_report.html` — combined tumour+normal QC report
- `results/alignment/COLO829/COLO829.dedup.bam` (+ `.bai`, + `.dedup.metrics.txt`)
- `results/alignment/COLO829BL/COLO829BL.dedup.bam` (+ `.bai`, + `.dedup.metrics.txt`)

Per the plan's Phase 1 verification step ("verify duplicate rates and coverage depth are appropriate for WGS somatic calling"): check the `.dedup.metrics.txt` `PERCENT_DUPLICATION` field and eyeball coverage via `samtools depth` or the MultiQC report before moving to Phase 2 — a pipeline that runs to completion isn't the same as one whose output looks right.
